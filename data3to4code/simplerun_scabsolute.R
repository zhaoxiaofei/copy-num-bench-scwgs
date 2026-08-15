#!/usr/bin/env Rscript
# simplerun_scabsolute.R <bam_dir> <out.tsv> <scAbsolute_dir>
#                        [binSize=500] [genome=hg19] [threads] [cellcycle] [splitPerChromosome]
#
# Runs scAbsolute (Schneider et al., Genome Biol 2024;25:62) over every BAM in <bam_dir> and
# writes ploidy calls in the format read by ploidy_tools.py:
#     cell  ploidy  rpc  used_reads  failure_reason
# Work files sit next to the TSV: .readCounts.rds (binned cache), .rds (per-cell objects),
# .cells/ (resume cache), .hmm_work/ (TF checkpoints).
#
# Environment facts this wrapper exists to encode -- each one cost a whole run to find:
#  1 scAbsolute is sourced, not installed, and reads BASEDIR / species / genome as FREE variables
#    in globalenv(). Omit `species` and R finds GenomeInfoDb::species instead, then dies at the
#    very end with "comparison (==) is possible only for atomic and list types".
#  2 upstream R/load_dependencies.R is unused: it calls future::plan("multiprocess"), defunct.
#  3 segment() applies PELT through apply(), so it is single-threaded whatever the plan is.
#    Parallelism therefore has to be over cells: bin once, then one cell per worker. batchSize=1
#    is upstream's default and the bin filter comes from the annotations, so a cell's call does
#    not depend on which other cells are in the run.
#  4 computeScale() scales in python via reticulate (numpy/pandas/tensorflow), on the method=
#    "error" path too. Pin the interpreter and import UP FRONT -- a missing pandas otherwise
#    surfaces 43 min in -- and cap threads, or one TF per worker oversubscribes the machine.
#  5 R/segmentation.py targets tensorflow_probability <= 0.14, whose DeferredModule took an
#    `args_fn` transform. Since 0.15 that keyword is swallowed by **kwargs and every cell dies in
#    the step-7 HMM with "experimental_from_mean_dispersion() missing 2 required positional
#    arguments" -- after segmentation and scaling are already paid for. check_tfp_api() turns
#    that into a one-second startup failure. Real fix: scabsolute-tfp015-deferredmodule.patch.
#  6 segmentation.py sets TF thread counts at MODULE level, which TF forbids once a context
#    exists (computeScale creates one). Sourcing it during setup, before any context, sets them
#    while it is still legal; the later repeat assigns the same values, which TF short-circuits.
#  7 hmm_path must be an explicit writable directory, or upstream builds /<cell>.hmm and
#    manager.save() fails with PermissionDeniedError.
#
# Speed: setup (16 library() calls, 7 sys.source() files, the TF warm-up) is memoised per R
# process instead of repeated per cell; each worker is sent only its own slice of the binned
# object rather than the whole matrix; binning and per-cell results are cached so a re-run
# resumes. splitPerChromosome=yes is the one flag that moves wall clock by an order of
# magnitude (PELT is superlinear in bin count) -- it is off by default because it changes where
# segment boundaries may fall, i.e. the method, not just its implementation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("Usage: Rscript simplerun_scabsolute.R <bam_dir> <out.tsv> <scAbsolute_dir> ",
         "[binSize] [genome] [threads] [cellcycle] [splitPerChromosome]")
}

opt  <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
flag <- function(i) tolower(opt(i, "no")) %in% c("yes", "true", "1")

input_dir   <- args[1]
output_tsv  <- args[2]
BASEDIR     <- normalizePath(args[3])                    # global, note 1
binSize     <- as.numeric(opt(4, 500))                   # kb
genome      <- opt(5, "hg19")                            # global, note 1
species     <- "Human"                                   # global, note 1
nthreads    <- suppressWarnings(as.integer(opt(6, NA)))
doCellcycle <- flag(7)                                   # S-phase tests a ploidy benchmark drops
splitChr    <- flag(8)                                   # see the speed note above

PKGS <- c("reticulate", "QDNAseq", "Biobase", "BiocGenerics", "GenomicRanges",
          "Rsamtools", "dplyr", "readr", "digest", "IRanges", "MASS", "robustbase", "S4Vectors",
          "matrixStats", "ggplot2")
SRC  <- c("data/changepoint/wrap_PELT.R", "R/scSegment.R", "R/scAbsolute.R", "R/core.R",
          "R/mean-variance.R", "R/visualization.R", "R/cellcycle.R")
PY   <- c("numpy", "pandas", "tensorflow", "tensorflow_probability")   # note 4


## Helpers --------------------------------------------------------------------------------------

# The TSV is tab-separated and newline-terminated, so a multi-line failure_reason would corrupt
# every field after it, not just its own. Tracebacks go to the log instead.
one_line <- function(x) if (length(x) == 1L && !is.na(x)) trimws(gsub("[[:space:]]+", " ", x)) else x

stage <- function(label, expr) {
    t0 <- Sys.time(); value <- force(expr)
    message(sprintf("[scAbsolute] %-22s %7.1f min", label,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    value
}

# Create, then prove writability with a real probe: permission/mount/ACL problems are worth
# catching now rather than at the first checkpoint write. Returns an absolute path, so workers do
# not depend on their own working directory.
writable_dir <- function(path, label) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) stop("cannot create ", label, ": ", path, call. = FALSE)
    path  <- normalizePath(path, mustWork = TRUE)
    probe <- tempfile("wtest-", tmpdir = path)
    if (!isTRUE(tryCatch(suppressWarnings(file.create(probe)), error = function(e) FALSE)))
        stop(label, " is not writable: ", path, call. = FALSE)
    unlink(probe, force = TRUE)
    path
}

# Everything a process needs before it can call scAbsolute. Memoised through options(), which
# persist for the life of an R process but are not inherited by multisession workers -- so this
# runs once per worker, not once per cell, and is still correct under plan(sequential).
setup <- function(cfg) {
    if (!isTRUE(getOption("scabs.loaded"))) {
        if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {          # note 4: pin the interpreter
            cand <- c(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python"),
                      Sys.which("python3"), Sys.which("python"))
            cand <- cand[nzchar(cand) & file.exists(cand)]
            if (length(cand)) Sys.setenv(RETICULATE_PYTHON = cand[[1]])
        }
        n <- as.character(cfg$cores)
        do.call(Sys.setenv, c(
            stats::setNames(as.list(rep(n, 6L)),
                            c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                              "NUMEXPR_NUM_THREADS", "TF_NUM_INTRAOP_THREADS",
                              "TF_NUM_INTEROP_THREADS")),
            list(TF_CPP_MIN_LOG_LEVEL = "2")))
        if (identical(Sys.getenv("scAbsoluteUseGPU"), "1")) {
            Sys.setenv(TF_FORCE_GPU_ALLOW_GROWTH = "true")       # several workers, one device
        } else {
            Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")              # the SVI fit is small; CPU is safer
        }

        suppressPackageStartupMessages(for (p in PKGS) library(p, character.only = TRUE))
        for (v in c("BASEDIR", "species", "genome")) assign(v, cfg[[v]], envir = globalenv())
        for (f in SRC) sys.source(file.path(cfg$BASEDIR, f), envir = globalenv())
        options(scabs.loaded = TRUE)
    }

    # Note 6. Fatal on failure: continuing only postpones the same problem to the per-cell HMM.
    if (isTRUE(cfg$tf) && !isTRUE(getOption("scabs.tf"))) {
        tryCatch(reticulate::source_python(file.path(cfg$BASEDIR, "R", "segmentation.py"),
                                           convert = TRUE),
                 error = function(e) stop("tensorflow warm-up failed: ", conditionMessage(e),
                                          call. = FALSE))
        options(scabs.tf = TRUE)
    }
    invisible(TRUE)
}

fail <- function(why) structure(list(message = why, traceback = NULL), class = "scAbsoluteFailure")

# Every per-cell value must be a non-NULL object of exactly one cell. scAbsolute returns NULL for
# a cell it drops in its own QC, and future.apply's gather then FLATTENS THAT AWAY: 12 cells go
# in, 11 values come back, and the whole run dies at the merge with
#     "the total number of elements (= 11) does not match ... in 'X' (= 12)"
# after every other cell has already been paid for. One silently dropped cell must cost that
# cell, not the run, so anything unusable becomes an ordinary scAbsoluteFailure here.
as_result <- function(res, why) {
    if (inherits(res, "scAbsoluteFailure")) return(res)
    if (is.null(res)) return(fail(why))
    n <- tryCatch(ncol(res), error = function(e) NA_integer_)
    if (length(n) == 1L && !is.na(n) && n < 1L) return(fail(paste0(why, " (0 cells in the result)")))
    res
}

# conditionMessage() on a reticulate error gives only the last line of the python exception --
# what went wrong, never where, and "where" is the whole question several frames into scAbsolute's
# own python. The traceback is per-process and gone once the worker exits, so read it in there.
py_traceback <- function() {
    e <- tryCatch(reticulate::py_last_error(), error = function(...) NULL)
    if (is.null(e)) return(NULL)
    # Newer reticulate returns an object whose print method formats the traceback; older versions
    # a plain list. Cover both.
    txt <- tryCatch(paste(utils::capture.output(print(e)), collapse = "\n"), error = function(...) "")
    if (!nzchar(trimws(txt)))
        txt <- paste(unlist(e[intersect(c("type", "value", "message", "traceback"), names(e))]),
                     collapse = "\n")
    if (nzchar(trimws(txt))) txt else NULL
}

# Import rather than merely look for: pip reports site-packages, scAbsolute needs `import` to
# succeed inside this R process, and a libstdc++/BLAS clash breaks the second without the first.
check_python <- function() {
    cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
    err <- vapply(PY, function(m) tryCatch({ reticulate::import(m); "" },
                                           error = function(e) conditionMessage(e)), character(1))
    message("[scAbsolute] python: ", if (is.null(cfg)) "<none found>" else cfg$python)
    bad <- names(err)[nzchar(err)]
    if (length(bad)) stop(sprintf(paste0(
        "python module(s) unusable inside R: %s\n  interpreter: %s\n",
        "  ModuleNotFoundError below => conda run -n scabsolute pip install %s\n",
        "  any other error => installed but unimportable in R; compare with\n",
        "      conda run -n scabsolute python -c 'import %s'\n%s"),
        paste(bad, collapse = ", "), if (is.null(cfg)) "?" else cfg$python,
        paste(bad, collapse = " "), bad[[1]],
        paste(sprintf("  --- import %s failed ---\n%s", bad, sub("\n+$", "", err[bad])),
              collapse = "\n")), call. = FALSE)
    invisible(TRUE)
}

# Note 5. Deliberately narrow: fires only when segmentation.py still passes args_fn AND the
# installed DeferredModule no longer accepts it. Anything undeterminable counts as fine, because
# refusing a good run is worse than reporting a bad one late. Introspection, not a version
# string: a patched checkout, a backport and a fork all read the same either way.
check_tfp_api <- function(basedir) {
    seg <- file.path(basedir, "R", "segmentation.py")
    if (!file.exists(seg) || !any(grepl("args_fn[[:space:]]*=", readLines(seg, warn = FALSE))))
        return(invisible(TRUE))
    par <- tryCatch(names(reticulate::import("inspect")$signature(
        reticulate::import("tensorflow_probability")$experimental$util$DeferredModule)$parameters),
        error = function(e) NULL)
    # "kwargs" is the positive control: both signatures end in *args/**kwargs, so a result
    # without it is not a parameter list and cannot answer the question.
    if (is.null(par) || !("kwargs" %in% par) || "args_fn" %in% par) return(invisible(TRUE))
    stop(sprintf(paste0(
        "scAbsolute's R/segmentation.py needs DeferredModule's `args_fn`, removed in ",
        "tensorflow_probability 0.15.0.\n  installed tfp: %s\n  DeferredModule now takes: (%s)\n",
        "  Left alone, every cell dies in the step-7 HMM AFTER segmentation and scaling, with\n",
        "      TypeError: experimental_from_mean_dispersion() missing 2 required positional arguments\n",
        "  Fix by patching %s (preferred: the rest of the environment works), or pin the old\n",
        "  API:  conda run -n scabsolute pip install 'tensorflow_probability<0.15'"),
        tryCatch(as.character(reticulate::import("tensorflow_probability")$`__version__`),
                 error = function(e) "?"),
        paste(par, collapse = ", "), seg), call. = FALSE)
}


## Preflight ------------------------------------------------------------------------------------

bams <- sort(Sys.glob(file.path(input_dir, "*.bam")))
if (!length(bams)) stop("No BAM file was found in ", input_dir)

ncores   <- as.integer(future::availableCores())
if (is.na(nthreads) || nthreads < 1L) nthreads <- max(1L, min(length(bams), ncores))

out_dir    <- writable_dir(dirname(output_tsv), "output directory")
output_tsv <- file.path(out_dir, basename(output_tsv))          # absolute, for the workers
hmm_root   <- writable_dir(paste0(output_tsv, ".hmm_work"), "HMM work directory")   # note 7
cell_dir   <- writable_dir(paste0(output_tsv, ".cells"), "per-cell cache directory")

# Allow the per-part QDNAseq objects to be exported to multisession workers.
options(future.globals.maxSize = 8 * 1024^3)

# parallelly refuses more workers than availableCores(), which under-reports on shared or
# containerised hosts (it follows cgroup quotas and the scheduler's variables). An explicit
# request is honoured by lifting the cap rather than by silently running narrower.
if (nthreads > ncores) {
    message(sprintf("[scAbsolute] availableCores() reports %d, %d requested: raising the limit",
                    ncores, nthreads))
    options(parallelly.maxWorkers.localhost = ceiling(nthreads / max(1L, ncores)) + 1)
}
if (nthreads > 1L) future::plan(future::multisession, workers = nthreads) else future::plan(future::sequential)

cfg <- list(BASEDIR = BASEDIR, species = species, genome = genome, tf = TRUE,
            cores = max(1L, ncores %/% max(1L, nthreads)))

message(sprintf(paste0("[scAbsolute] %d cell(s), binSize = %g kb, genome = %s, %d worker(s) of ",
                       "%d core(s), cellcycle = %s, splitPerChromosome = %s"),
                length(bams), binSize, genome, nthreads, ncores, doCellcycle, splitChr))

setup(modifyList(cfg, list(tf = FALSE)))     # the parent bins; only workers need TensorFlow
check_python()
check_tfp_api(BASEDIR)
if (!("hmm_path" %in% names(formals(scAbsolute))))
    stop("the loaded scAbsolute() has no `hmm_path` argument; refusing to run, see note 7",
         call. = FALSE)


## 1. Bin once, cached --------------------------------------------------------------------------

cache_file <- paste0(output_tsv, ".readCounts.rds")
key        <- digest::digest(list(basename(bams), file.info(bams)$size, binSize, genome, species))
readCounts <- NULL

if (file.exists(cache_file)) {
    hit <- try(readRDS(cache_file), silent = TRUE)
    if (!inherits(hit, "try-error") && identical(hit$key, key)) {
        readCounts <- hit$readCounts
        message("[scAbsolute] reusing the binned read counts in ", cache_file)
    }
    rm(hit)
}

if (is.null(readCounts)) {
    readCounts <- stage("readData", readData(bams, binSize = binSize, species = species,
                                             genome = genome, filterChromosomes = c("MT")))
    saveRDS(list(key = key, readCounts = readCounts), cache_file)
}

cells <- sub("\\.bam$", "", as.character(Biobase::pData(readCounts)[["name"]]))
stopifnot(length(cells) == length(bams), !anyNA(cells))   # never write a misaligned TSV


## 2. Scale each cell to absolute copy number ----------------------------------------------------

# Cells are dealt round-robin into chunks and each chunk is its own future, so slow cells spread
# across workers and finished workers pick up the queue. Each chunk carries only its own slice of
# the binned object: the full matrix is never serialised to every worker.
ix     <- seq_along(cells)
groups <- unname(split(ix, rep(seq_len(min(length(ix), max(1L, 4L * nthreads))), length.out = length(ix))))
parts  <- lapply(groups, function(g) list(ix = g, cells = cells[g], rc = readCounts[, g, drop = FALSE]))
rm(readCounts); invisible(gc(FALSE))

tag <- substr(key, 1L, 8L)

run_part <- function(part) {
    ready <- FALSE                            # lazy: an already-cached chunk never pays setup
    lapply(seq_along(part$cells), function(j) {
        cell <- part$cells[j]
        done <- file.path(cell_dir, paste0(tag, "_", cell, ".rds"))

        # Segmentation is the expensive step; a downstream failure must not buy it twice. Delete
        # the .cells directory to force a clean run.
        if (file.exists(done)) {
            hit <- try(readRDS(done), silent = TRUE)
            if (!inherits(hit, "try-error")) {
                # An earlier run could have cached a NULL, which would resurrect the gather
                # failure on every future run; drop it and recompute once instead.
                if (!inherits(as_result(hit, ""), "scAbsoluteFailure")) return(hit)
                unlink(done)
            }
        }

        if (!ready) { setup(cfg); ready <<- TRUE }

        hmm <- file.path(hmm_root, cell)      # note 7; the root is already write-probed
        dir.create(hmm, recursive = TRUE, showWarnings = FALSE)
        set.seed(2020)                        # upstream's seed, per cell, so a cell's call does
                                              # not depend on its position in the run
        # get() rather than a bare name, so the globals scanner does not ship the parent's copy
        # of scAbsolute (and its dependencies) to every worker that already sourced them.
        sca <- get("scAbsolute", envir = globalenv())

        # A dropped cell raises nothing, so its warnings are the only evidence of why it went.
        # Recorded, not muffled: they still reach the log as usual.
        warns <- character(0)

        res <- withCallingHandlers(
            tryCatch(
                sca(part$rc[, j, drop = FALSE], binSize = binSize, species = species, genome = genome,
                    minPloidy = 1.5, maxPloidy = 6.0, ploidyWindow = 0.1,  # bound the search, not the answer
                    splitPerChromosome = splitChr, doCellcycle = doCellcycle,
                    quick = TRUE,             # drops the mean-variance model: descriptive columns
                    readPositionModel = FALSE, # only, nothing that feeds ploidy or rpc
                    hmm_path = hmm),          # CRITICAL: never NULL/empty
                error = function(e) structure(list(message = conditionMessage(e),
                                                   traceback = py_traceback()),
                                              class = "scAbsoluteFailure")),
            warning = function(w) warns <<- c(warns, conditionMessage(w)))

        res <- as_result(res, paste0(
            "scAbsolute returned no object for this cell",
            if (length(warns)) paste0("; warnings: ", paste(utils::tail(warns, 3L), collapse = " | ")) else ""))

        if (!inherits(res, "scAbsoluteFailure")) try(saveRDS(res, done), silent = TRUE)
        res
    })
}

# Use futures directly and normalise every worker return before merging.  A dropped cell may
# surface as NULL at this boundary; it is a per-cell failure, not a reason to abort the whole run.
chunks <- stage("scAbsolute (all cells)", {
    fs <- lapply(parts, function(part) future::future(run_part(part), seed = TRUE))
    Map(function(f, part) {
        x <- tryCatch(future::value(f),
                      error = function(e) fail(paste0("worker failed: ", conditionMessage(e))))
        if (inherits(x, "scAbsoluteFailure")) return(rep(list(x), length(part$cells)))
        if (!is.list(x) || length(x) != length(part$cells))
            return(rep(list(fail(sprintf("worker returned %d result(s) for %d cell(s)",
                                         length(x), length(part$cells)))), length(part$cells)))
        Map(function(y, cell) as_result(y, paste0("worker returned no object for ", cell)),
            x, part$cells)
    }, fs, parts)
})

results <- vector("list", length(cells))      # back into BAM order, whatever the chunking was
for (k in seq_along(parts)) results[parts[[k]]$ix] <- chunks[[k]]
results <- Map(function(x, cell) as_result(x, paste0("no result collected for ", cell)),
               results, cells)


## 3. Collect the per-cell ploidy calls ----------------------------------------------------------

# Older scAbsolute versions do not record every field, so each one is read defensively.
field <- function(pd, name, cast) if (name %in% colnames(pd)) cast(pd[[name]][1]) else cast(NA)

row_of <- function(res, cell) {
    res <- as_result(res, paste0("no result available for ", cell))
    if (inherits(res, "scAbsoluteFailure"))
        return(data.frame(cell = cell, ploidy = NA_real_, rpc = NA_real_, used_reads = NA_real_,
                          failure_reason = one_line(paste0("error: ", res$message)),
                          stringsAsFactors = FALSE))
    pd <- Biobase::pData(res)
    nm <- field(pd, "name", as.character)
    data.frame(cell           = if (is.na(nm)) cell else sub("\\.bam$", "", nm),
               ploidy         = field(pd, "ploidy", as.numeric),
               rpc            = field(pd, "rpc", as.numeric),
               used_reads     = field(pd, "used.reads", as.numeric),
               failure_reason = one_line(field(pd, "failure_reason", as.character)),
               stringsAsFactors = FALSE)
}

out <- do.call(rbind, Map(row_of, results, cells))

# A cell that failed to scale carries a meaningless ploidy: blank it so it counts as unusable
# downstream rather than scoring as a wrong answer. "" is not a failure reason -- upstream writes
# NA_character_ for a cell that passed, and "" would otherwise condemn the whole run.
failed <- (!is.na(out$failure_reason) & nzchar(out$failure_reason)) | !is.finite(out$rpc) | out$rpc <= 0
out$ploidy[failed] <- NA

# One line per failure in the TSV, the whole traceback in the log: the last line names the
# symptom, the frames underneath it name the call site.
for (i in which(failed)) {
    message("[scAbsolute] FAILED: ", cells[i],
            if (!is.na(out$failure_reason[i])) paste0(": ", out$failure_reason[i]) else "")
    tb <- if (inherits(results[[i]], "scAbsoluteFailure")) results[[i]]$traceback else NULL
    if (!is.null(tb)) message(paste0("  | ", strsplit(tb, "\n", fixed = TRUE)[[1]], collapse = "\n"))
}

ok <- !vapply(results, inherits, logical(1), "scAbsoluteFailure")
if (any(ok)) saveRDS(results[ok], paste0(output_tsv, ".rds"))

write.table(out, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message(sprintf("[scAbsolute] ploidy calls written to %s (%d of %d cells succeeded)",
                output_tsv, sum(!failed), nrow(out)))

if (all(failed)) quit(status = 1)
