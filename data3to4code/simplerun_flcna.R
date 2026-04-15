#!/usr/bin/env Rscript
# simplerun_flcna.R
# Usage:
#   Rscript simplerun_flcna.R <input_bam_dir> <ref_fasta> <output_csv> <mappability_bigwig>
#
# Pipeline:
#   1. Build 200 kb bins from FASTA .fai
#   2. Count reads in BAMs with bedtools multicov
#   3. Compute GC with bedtools nuc
#   4. Compute mappability from bigWig on the SAME bins
#   5. Run FLCNA

# https://www.doubao.com/chat/38418713917526274 : bin-size of 100kb requires too much (200GB+) RAM
bin_size <- 500*1000 # NOTE: if 200kb takes too much runtime and/or RAM, set it to 500kb million.

# Set console width to a very large value (eliminates line wrapping/truncation)
options(width = 10000)

# Max out warning message length (R's hard maximum is 8172 characters)
options(warning.length = 5000)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript simplerun_flcna.R <input_bam_dir> <ref_fasta> <output_csv> <mappability_bigwig>")
}

input_dir      <- args[1]
ref_fasta      <- args[2]
output_csv     <- args[3]
wgECMA_bigwig  <- args[4]

# https://www.doubao.com/chat/38421245670502402 
# how to prevent soft exit (with exit_code=0) from recursive errors

suppressPackageStartupMessages({
    library(FLCNA)
    library(GenomicRanges)
    library(IRanges)
    library(rtracklayer)
})

# <FIX_MISSING_CNV_STATE, DEBUG_FROM='https://www.doubao.com/chat/38421178395960578'>

# Step 1: Capture original FLCNA function
original_Para_init <- FLCNA:::Para_init

# Step 2: Wrapper patch - Linear INTERPOLATION + EXTRAPOLATION (VALID STATES UNMODIFIED)
modified_Para_init <- function(...) {
  # Run original function (raw output)
  out <- original_Para_init(...)
  n_states <- length(out$p) # Total CNV states (1-5)
  state_order <- 1:n_states  # CNV state positions (1=lowest, 5=highest)

  # 1. Identify VALID states (NO CHANGES to these values)
  valid_idx <- which(
    !is.na(out$mu) & !is.infinite(out$mu) &
    !is.na(out$sigma) & out$sigma > 0 &
    out$p > 0 & !is.na(out$p)
  )
  invalid_idx <- setdiff(state_order, valid_idx)

  if (length(invalid_idx) > 0) {
    # Extract sorted valid data (preserve original valid mu/sigma/p)
    valid_x <- sort(valid_idx)          # Valid CNV state positions
    valid_y <- out$mu[valid_x]          # Valid log2 ratios (UNCHANGED)

    # 2. LINEAR INTERPOLATION + EXTRAPOLATION for invalid mu
    # Uses valid states to fill gaps (interpolate) and edges (extrapolate)
    out$mu[invalid_idx] <- approx(
      x = valid_x,
      y = valid_y,
      xout = invalid_idx,
      method = "linear",
      rule = 2  # rule=2 enables EXTAPOLATION (matches valid trend)
    )$y

    # 3. Fill invalid sigma/p (VALID values stay 100% original)
    out$sigma[invalid_idx] <- mean(out$sigma[valid_idx], na.rm = TRUE)
    out$p[invalid_idx] <- 0.001
  }

  # Final cleanup (no changes to valid values)
  out$sigma <- pmax(out$sigma, 1e-6)
  out$p <- pmax(out$p, 1e-6)
  out$p <- out$p / sum(out$p) # Normalize probabilities

  # Enforce: Log2ratio INCREASES with CNV state (1→5)
  reorder <- order(out$mu)
  out$mu <- out$mu[reorder]
  out$sigma <- out$sigma[reorder]
  out$p <- out$p[reorder]

  return(out)
}

# Step 3: Apply the monkey-patch
assignInNamespace("Para_init", modified_Para_init, ns = "FLCNA")

# </FIX_MISSING_CNV_STATE>

# -----------------------------
# 1. Find BAM files
# -----------------------------
bam_files <- sort(Sys.glob(file.path(input_dir, "*.bam")))
bam_files <- bam_files[!grepl("\\.bam\\.bai$", bam_files)]

if (length(bam_files) == 0) {
    stop("No BAM files found in ", input_dir)
}

cat("Found", length(bam_files), "BAM files\n")

# -----------------------------
# 2. Create bins from .fai
# -----------------------------
fai_file <- paste0(ref_fasta, ".fai")
if (!file.exists(fai_file)) {
    stop("FAI index not found: ", fai_file)
}

fai <- read.table(
    fai_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)
colnames(fai) <- c("chr", "length", "offset", "bases_per_line", "bytes_per_line")

# https://gemini.google.com/app/e9816bb5aa8a4328 : FLCNA stops with runtime error if incorporating non-autosomes
allowed_chrs <- paste0("chr", c(1:22))
fai <- fai[fai$chr %in% allowed_chrs, , drop = FALSE]

if (nrow(fai) == 0) {
    stop("No allowed chromosomes found in FAI. Expected names like chr1..chr22, chrX, chrY")
}

bed_file <- paste0(output_csv, ".bed.tmp")

bed_df_list <- vector("list", nrow(fai))
for (i in seq_len(nrow(fai))) {
    chr_name <- fai$chr[i]
    chr_len  <- fai$length[i]

    starts <- seq(0, chr_len - 1, by = bin_size)
    ends   <- pmin(starts + bin_size, chr_len)

    bed_df_list[[i]] <- data.frame(
        chr   = chr_name,
        start = starts,
        end   = ends,
        stringsAsFactors = FALSE
    )
}

bed_df <- do.call(rbind, bed_df_list)
bed_df[, 2] <- as.integer(bed_df[, 2])
bed_df[, 3] <- as.integer(bed_df[, 3])
write.table(
    bed_df,
    file = bed_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

cat("Generated", nrow(bed_df), "bins\n")

# -----------------------------
# 3. Generate read counts
# -----------------------------
multicov_out <- paste0(output_csv, ".cov.tmp")

cmd_multicov <- c(
    "multicov",
    "-bams", bam_files,
    "-bed", bed_file
    # "-o", multicov_out
)

cat("Running: bedtools", paste(cmd_multicov, collapse = " "), "\n")

# --- FIX 2: Run bedtools and check exit status ---
multicov_status <- system2(
    command = "bedtools",
    args = cmd_multicov,
    stdout = multicov_out,
    stderr = paste0(multicov_out, '.stderr') # TRUE
)

# Check if bedtools itself failed (exit code != 0)
exit_status <- attr(multicov_status, "status")
if (!is.null(exit_status) && exit_status != 0) {
    stop(
        "bedtools multicov failed (exit status ", exit_status, ")\n",
        "bedtools stderr:\n", paste(multicov_status, collapse = "\n")
    )
}

# --- FIX 3: Check output file ---
if (!file.exists(multicov_out)) {
    stop("bedtools did not create output file: ", multicov_out)
}
if (file.info(multicov_out)$size == 0) {
    stop(
        "bedtools created an empty file. Common reasons:\n",
        "1. BAM chromosomes (e.g., '1') don't match BED (e.g., 'chr1')\n",
        "2. BAM files have zero mapped reads (check with 'samtools flagstat')"
    )
}

if (!file.exists(multicov_out) || file.info(multicov_out)$size == 0) {
    stop("bedtools multicov failed or produced empty output: ", multicov_out)
}

multicov <- read.table(
    multicov_out,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)

n_bams <- length(bam_files)
expected_cols <- 3 + n_bams
if (ncol(multicov) != expected_cols) {
    stop(
        "Unexpected multicov output: got ", ncol(multicov),
        " columns, expected ", expected_cols,
        " (3 BED columns + ", n_bams, " BAM count columns)"
    )
}

Y_raw <- as.matrix(multicov[, 4:(3 + n_bams), drop = FALSE])
storage.mode(Y_raw) <- "numeric"

colnames(Y_raw) <- sapply(
    bam_files,
    function(f) tools::file_path_sans_ext(basename(f))
)
rownames(Y_raw) <- paste0(multicov[, 1], ":", multicov[, 2], "-", multicov[, 3])

cat("Y_raw dimensions:", nrow(Y_raw), "x", ncol(Y_raw), "\n")

# -----------------------------
# 4. Generate GC content using bedtools nuc
# -----------------------------
gc_out <- paste0(output_csv, ".gc.tmp")

cmd_gc <- c(
    "nuc",
    "-fi", ref_fasta,
    "-bed", bed_file
    # "-o", gc_out
)

cat("Running: bedtools", paste(cmd_gc, collapse = " "), "\n")

if (TRUE) {
gc_status <- system2(
    command = "bedtools",
    args = cmd_gc,
    stdout = gc_out,
    stderr = paste0(gc_out, '.stderr') # TRUE
)
}

if (!file.exists(gc_out) || file.info(gc_out)$size == 0) {
    stop("bedtools nuc failed or produced empty output: ", gc_out)
}

gc_data <- read.table(
    gc_out,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = ""
)

gc_col <- grep("pct_gc", colnames(gc_data), value = TRUE)
if (length(gc_col) == 0) {
    stop("Could not find pct_gc column in bedtools nuc output")
}

gc_pct <- gc_data[[gc_col[1]]] * 100

if (length(gc_pct) != nrow(multicov)) {
    stop(
        "GC vector length mismatch: length(gc_pct) = ", length(gc_pct),
        ", nrow(multicov) = ", nrow(multicov)
    )
}

# -----------------------------
# 5. Build ref_raw from SAME bins as multicov
# -----------------------------
# Convert BED-style coordinates to GRanges:
# BED: 0-based, half-open
# GRanges: 1-based, closed
bins <- GRanges(
    seqnames = multicov[, 1],
    ranges = IRanges(
        start = multicov[, 2] + 1,
        end   = multicov[, 3]
    )
)

# Preserve seqlevels order
seqlevels(bins) <- unique(as.character(seqnames(bins)))

# Attach GC from bedtools nuc computed on the same bins
bins$gc <- as.numeric(gc_pct)

# Compute mean mappability on the same bins
cat("Loading mappability bigWig:", wgECMA_bigwig, "\n")
if (!file.exists(wgECMA_bigwig) && !grepl("^https?://", wgECMA_bigwig)) {
    stop("Mappability bigWig not found: ", wgECMA_bigwig)
}

bw_file <- BigWigFile(wgECMA_bigwig)
mapp_summary <- summary(bw_file, which = bins, type = "mean")

if (FALSE) {
# rtracklayer::summary may return a DataFrame-like object with score column
if ("score" %in% colnames(as.data.frame(mapp_summary))) {
    bins$mapp <- as.numeric(mapp_summary$score)
} else {
    bins$mapp <- as.numeric(mapp_summary)
}
}

# <FIX>
cat("Loading mappability bigWig:", wgECMA_bigwig, "\n")
if (!file.exists(wgECMA_bigwig) && !grepl("^https?://", wgECMA_bigwig)) {
    stop("Mappability bigWig not found: ", wgECMA_bigwig)
}

bw_file <- BigWigFile(wgECMA_bigwig)

# Import mappability over the exact bins
mapp_list <- import(bw_file, which = bins, as = "NumericList")

# One mean mappability value per bin
bins$mapp <- vapply(
    mapp_list,
    function(x) {
        if (length(x) == 0 || all(is.na(x))) {
            NA_real_
        } else {
            mean(x, na.rm = TRUE)
        }
    },
    numeric(1)
)

# Replace missing values
bins$mapp[is.na(bins$mapp)] <- 0

cat("length(bins$mapp):", length(bins$mapp), "\n")
stopifnot(length(bins$mapp) == length(bins))

# </FIX>

# Replace missing mappability values
bins$mapp[is.na(bins$mapp)] <- 0

ref_raw <- bins

cat("length(ref_raw):", length(ref_raw), "\n")
if (nrow(Y_raw) != length(ref_raw)) {
    stop(
        "Internal mismatch before FLCNA_QC: nrow(Y_raw) = ", nrow(Y_raw),
        ", length(ref_raw) = ", length(ref_raw)
    )
}

# -----------------------------
# 6. Run FLCNA pipeline
# -----------------------------
cat("Running FLCNA QC...\n")
QCobject <- FLCNA_QC(
    Y_raw = Y_raw,
    ref_raw = ref_raw,
    mapp_thresh = 0.9,
    gc_thresh = c(20, 80)
)

cat("Running FLCNA normalization...\n")
log2Rdata <- FLCNA_normalization(
    Y   = QCobject$Y,
    gc  = QCobject$ref$gc,
    map = QCobject$ref$mapp
)

cat("Running FLCNA CNA detection...\n")
K_values <- if (n_bams >= 6) c(2, 3, 4) else c(2, 3)

output_FLCNA <- FLCNA(
    K = K_values,
    lambda = 3,
    y = t(log2Rdata),
    ref = QCobject$ref
)

cat("Saving workspace for potential debugging...\n")
save(output_FLCNA, log2Rdata, QCobject, file = paste0(output_csv, ".pre_cna.RData"))

cat("Running FLCNA CNA clustering...\n")
# Wrap in tryCatch to ensure we get a clean error if it fails
CNA_output <- tryCatch({
    CNA.out(
        mean.matrix = output_FLCNA$mu.hat.best,
        Clusters    = output_FLCNA$s.hat.best,
        LRR         = log2Rdata,
        QC_ref      = QCobject$ref,
        cutoff      = 0.80,
        L           = 100
    )
}, error = function(e) {
    message("CNA.out failed! You can debug this by loading: ", paste0(output_csv, ".pre_cna.RData"))
    print(e)
    q(save = "no", status = 1)  # Force R to exit with non-zero status
})

# -----------------------------
# 7. Write output
# -----------------------------
write.csv(CNA_output, file = output_csv, row.names = FALSE)
cat("FLCNA output written to:", output_csv, "\n")

# -----------------------------
# 8. Optional cleanup
# -----------------------------
# unlink(c(bed_file, multicov_out, gc_out))

