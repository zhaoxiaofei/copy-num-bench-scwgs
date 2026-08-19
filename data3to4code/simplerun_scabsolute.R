#!/usr/bin/env Rscript
# Usage:
#   Rscript simplerun_scabsolute.R <bam_dir> <out.tsv> <scAbsolute_dir> \
#       [binSize=500] [genome=hg19] [threads]
#
# Examples:
#   Rscript simplerun_scabsolute.R bams calls.tsv scAbsolute 500 hg19 16
#   SCABS_THREADS=16 Rscript simplerun_scabsolute.R bams calls.tsv scAbsolute 500 hg19
#
# Parallelism is across BAMs. Each BAM worker is CPU-only and pinned to one
# logical CPU to avoid TensorFlow/BLAS oversubscription.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript simplerun_scabsolute.R",
      "<bam_dir> <out.tsv> <scAbsolute_dir> [binSize] [genome] [threads]"
    ),
    call. = FALSE
  )
}

bam_dir <- args[[1]]
out_tsv <- args[[2]]
BASEDIR <- normalizePath(args[[3]], mustWork = TRUE)
binSize <- if (length(args) >= 4L) as.numeric(args[[4]]) else 500
genome  <- if (length(args) >= 5L) args[[5]] else "hg19"
species <- if (genome %in% c("mm10", "GRCm38")) "Mouse" else "Human"

if (!is.finite(binSize) || binSize <= 0) {
  stop("binSize must be a positive number", call. = FALSE)
}

# IMPORTANT: this must happen before the first Python/TensorFlow import.
#
# scAbsolute's computeScale() imports data/scAbsolute.py before the later
# segmentation.py is imported. segmentation.py itself tries to set
# CUDA_VISIBLE_DEVICES=-1, but that is too late if TensorFlow was already
# initialized by computeScale(). Running BAMs in parallel otherwise makes all
# children compete for GPU 0 and can kill workers with CUDA OOM.
Sys.setenv(
  CUDA_VISIBLE_DEVICES = "-1",
  TF_CPP_MIN_LOG_LEVEL = "2",

  # These constrain BLAS/OpenMP and are also useful TensorFlow hints.
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  TF_NUM_INTRAOP_THREADS = "1",
  TF_NUM_INTEROP_THREADS = "1"
)

bams <- sort(Sys.glob(file.path(bam_dir, "*.bam")))
if (!length(bams)) {
  stop("No BAM files found in ", bam_dir, call. = FALSE)
}

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(
  normalizePath(dirname(out_tsv), mustWork = TRUE),
  basename(out_tsv)
)

# Load the same source-based scAbsolute code as the original runner.
pkgs <- c(
  "reticulate", "QDNAseq", "Biobase", "BiocGenerics", "GenomicRanges",
  "Rsamtools", "dplyr", "readr", "IRanges", "MASS", "robustbase",
  "S4Vectors", "matrixStats", "ggplot2"
)

missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop(
    "Missing R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  invisible(lapply(pkgs, library, character.only = TRUE))
)

src <- c(
  "data/changepoint/wrap_PELT.R",
  "R/scSegment.R",
  "R/scAbsolute.R",
  "R/core.R",
  "R/mean-variance.R",
  "R/visualization.R",
  "R/cellcycle.R"
)

src_paths <- file.path(BASEDIR, src)
if (any(!file.exists(src_paths))) {
  stop(
    "Missing scAbsolute source file(s): ",
    paste(src[!file.exists(src_paths)], collapse = ", "),
    call. = FALSE
  )
}
invisible(lapply(src_paths, source))

# Number of concurrent BAMs:
# 1) explicit sixth argument
# 2) SCABS_THREADS
# 3) scheduler allocation
# 4) conservative default of <= 64
requested <- NA_integer_

if (length(args) >= 6L && nzchar(args[[6]])) {
  requested <- suppressWarnings(as.integer(args[[6]]))
}

if (is.na(requested) || requested < 1L) {
  requested <- suppressWarnings(
    as.integer(Sys.getenv("SCABS_THREADS", ""))
  )
}

if (is.na(requested) || requested < 1L) {
  scheduler <- suppressWarnings(as.integer(c(
    Sys.getenv("SLURM_CPUS_PER_TASK", ""),
    Sys.getenv("NSLOTS", ""),
    Sys.getenv("PBS_NP", "")
  )))
  scheduler <- scheduler[is.finite(scheduler) & scheduler > 0L]

  if (length(scheduler)) {
    requested <- scheduler[[1]]
  } else {
    detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (length(detected) != 1L || is.na(detected) || detected < 1L) {
      detected <- 1L
    }
    requested <- min(64L, detected)
  }
}

nworkers <- max(1L, min(length(bams), requested))

# CPU affinity is the hard limit on actual CPU use. TensorFlow's older v1
# session in data/scAbsolute.py uses n_cores=0, which means automatic thread
# selection. Environment variables alone therefore do not reliably guarantee
# one active core per child. On Linux, mcaffinity() restricts every thread in
# the child process to the selected CPU.
allowed_cpus <- tryCatch(parallel::mcaffinity(), error = function(e) NULL)
can_pin <- !is.null(allowed_cpus) && length(allowed_cpus) >= 1L

if (can_pin && nworkers > length(allowed_cpus)) {
  message(sprintf(
    "[scAbsolute] requested %d workers but only %d CPUs are in this process's affinity mask; using %d workers",
    nworkers, length(allowed_cpus), length(allowed_cpus)
  ))
  nworkers <- length(allowed_cpus)
}

message(sprintf(
  "[scAbsolute] %d BAM(s), %d parallel worker(s), GPU disabled, CPU pinning = %s",
  length(bams), nworkers, if (can_pin) "yes" else "no"
))

if (can_pin) {
  message(
    "[scAbsolute] allowed logical CPUs: ",
    paste(allowed_cpus, collapse = ",")
  )
}

# Each BAM gets its own HMM checkpoint directory.
hmm_root <- paste0(out_tsv, ".hmm_work")
dir.create(hmm_root, recursive = TRUE, showWarnings = FALSE)

get1 <- function(pd, name, default) {
  if (name %in% names(pd) && length(pd[[name]]) >= 1L) {
    pd[[name]][[1]]
  } else {
    default
  }
}

clean_error <- function(x) {
  trimws(gsub("[\t\r\n]+", " ", as.character(x)))
}

run_one <- function(i) {
  bam  <- bams[[i]]
  cell <- sub("\\.bam$", "", basename(bam), ignore.case = TRUE)

  # Pin this child to exactly one logical CPU. The mapping cycles only after
  # all CPUs in the inherited affinity mask have been used.
  if (can_pin) {
    cpu <- allowed_cpus[[ ((i - 1L) %% length(allowed_cpus)) + 1L ]]
    got <- tryCatch(
      parallel::mcaffinity(cpu),
      error = function(e) NULL
    )
    if (is.null(got)) {
      warning("Could not set CPU affinity for ", cell)
    }
  }

  # Repeat the environment settings inside the forked child. Most importantly,
  # CUDA must be hidden before scAbsolute causes Python/TensorFlow to load.
  Sys.setenv(
    CUDA_VISIBLE_DEVICES = "-1",
    TF_CPP_MIN_LOG_LEVEL = "2",
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    TF_NUM_INTRAOP_THREADS = "1",
    TF_NUM_INTEROP_THREADS = "1"
  )

  hmm <- file.path(hmm_root, sprintf("%06d", i))
  unlink(hmm, recursive = TRUE, force = TRUE)
  dir.create(hmm, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(hmm, recursive = TRUE, force = TRUE), add = TRUE)

  message(sprintf(
    "[scAbsolute] [%d/%d] %s%s",
    i, length(bams), cell,
    if (can_pin) paste0(" [cpu ", cpu, "]") else ""
  ))

  tryCatch({
    set.seed(2020)

    scaledCN <- scAbsolute(
      bam,
      binSize = binSize,
      species = species,
      genome = genome,
      minPloidy = 1.5,
      maxPloidy = 6.0,
      ploidyWindow = 0.1,
      batchSize = 1,
      doCellcycle = FALSE,
      quick = TRUE,
      readPositionModel = FALSE,
      hmm_path = hmm
    )

    pd <- Biobase::pData(scaledCN)
    if (nrow(pd) < 1L) {
      stop("scAbsolute returned no cell metadata")
    }

    row <- data.frame(
      cell = sub(
        "\\.bam$", "",
        as.character(get1(pd, "name", basename(bam))),
        ignore.case = TRUE
      ),
      ploidy = suppressWarnings(
        as.numeric(get1(pd, "ploidy", NA_real_))
      ),
      rpc = suppressWarnings(
        as.numeric(get1(pd, "rpc", NA_real_))
      ),
      used_reads = suppressWarnings(
        as.numeric(get1(pd, "used.reads", NA_real_))
      ),
      failure_reason = as.character(
        get1(pd, "failure_reason", NA_character_)
      ),
      stringsAsFactors = FALSE
    )

    failed <- (
      !is.na(row$failure_reason) & nzchar(row$failure_reason)
    ) | !is.finite(row$rpc) | row$rpc <= 0

    row$ploidy[failed] <- NA_real_
    row
  }, error = function(e) {
    structure(
      list(
        cell = cell,
        message = clean_error(conditionMessage(e))
      ),
      class = "scAbsoluteParallelError"
    )
  })
}

set.seed(2020)

# mc.preschedule=FALSE intentionally gives each BAM a fresh forked R process.
# That is important for reticulate/TensorFlow: a worker that has initialized
# TensorFlow for one cell is not reused for another cell.
if (nworkers == 1L) {
  rows <- lapply(seq_along(bams), run_one)
} else {
  rows <- parallel::mclapply(
    seq_along(bams),
    run_one,
    mc.cores = nworkers,
    mc.preschedule = FALSE,
    mc.set.seed = TRUE
  )
}

is_error <- vapply(
  rows,
  function(x) inherits(x, "scAbsoluteParallelError"),
  logical(1)
)
is_missing <- !vapply(rows, is.data.frame, logical(1))
bad <- is_error | is_missing

if (any(bad)) {
  msgs <- vapply(which(bad), function(i) {
    x <- rows[[i]]
    if (inherits(x, "scAbsoluteParallelError")) {
      sprintf("%s: %s", x$cell, x$message)
    } else {
      sprintf(
        "%s: worker exited without returning a result",
        basename(bams[[i]])
      )
    }
  }, character(1))

  stop(
    sprintf(
      "%d BAM(s) failed:\n%s",
      length(msgs),
      paste(msgs, collapse = "\n")
    ),
    call. = FALSE
  )
}

# mclapply returns values in input order, preserving the old sorted-BAM output.
out <- do.call(rbind, rows)

unlink(hmm_root, recursive = TRUE, force = TRUE)

write.table(
  out,
  out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

failed <- (
  !is.na(out$failure_reason) & nzchar(out$failure_reason)
) | !is.finite(out$rpc) | out$rpc <= 0

message(sprintf(
  "[scAbsolute] wrote %s: %d/%d cells succeeded",
  out_tsv,
  sum(!failed),
  nrow(out)
))
