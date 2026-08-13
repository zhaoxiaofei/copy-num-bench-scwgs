#!/usr/bin/env Rscript
# simplerun_scabsolute.R
# Usage: Rscript simplerun_scabsolute.R <input_bam_dir> <output_tsv> <scAbsolute_dir> [binSize] [genome]
#
# Runs scAbsolute (Schneider et al., Genome Biol 2024;25:62) on every BAM file in
# input_bam_dir and writes the per-cell ploidy calls to output_tsv in the ploidy-calls format
# read by ploidy_tools.py: a header line, then one row per cell with the columns
#   cell, ploidy, rpc, used_reads, failure_reason
# The full QDNAseq object is kept next to the TSV as <output_tsv>.rds, because it also holds
# the per-bin absolute copy numbers and the replication-status calls.
#
# NOTE scAbsolute is not an R package: it is loaded by sourcing its R/ files, and those files
# read their reference data through a global BASEDIR, which is why it is set here rather than
# passed as an argument. The upstream R/load_dependencies.R is deliberately not used, because
# it calls future::plan("multiprocess"), which is defunct in current versions of the future
# package; everything else it does is replayed below.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript simplerun_scabsolute.R <input_bam_dir> <output_tsv> <scAbsolute_dir> [binSize] [genome]")
}

input_dir  <- args[1]
output_tsv <- args[2]
BASEDIR    <- normalizePath(args[3])                                      # global, see NOTE above
binSize    <- if (length(args) >= 4) as.numeric(args[4]) else 500         # kb
genome     <- if (length(args) >= 5) args[5] else "hg19"

suppressPackageStartupMessages({
    library(reticulate);   library(future.apply); library(QDNAseq)
    library(Biobase);      library(BiocGenerics); library(GenomicRanges)
    library(Rsamtools);    library(dplyr);        library(readr)
    library(digest);       library(IRanges);      library(MASS)
    library(robustbase);   library(S4Vectors);    library(matrixStats)
    library(ggplot2)
})
try(future::plan(future::multisession), silent = TRUE)
for (f in c("data/changepoint/wrap_PELT.R", "R/scSegment.R", "R/scAbsolute.R", "R/core.R",
            "R/mean-variance.R", "R/visualization.R", "R/cellcycle.R")) {
    source(file.path(BASEDIR, f))
}

bamfiles <- sort(Sys.glob(file.path(input_dir, "*.bam")))
if (length(bamfiles) == 0) stop(paste("No BAM file was found in", input_dir))
cat("Running scAbsolute on", length(bamfiles), "cells with binSize =", binSize, "kb\n")

# minPloidy/maxPloidy are scAbsolute's own defaults: they bound the search, not the answer.
cn <- scAbsolute(bamfiles, binSize = binSize, genome = genome,
                 minPloidy = 1.2, maxPloidy = 10.0, ploidyWindow = 0.1)
saveRDS(cn, paste0(output_tsv, ".rds"))

pd <- Biobase::pData(cn)
# Older scAbsolute versions do not record every field, so each one is read defensively.
column <- function(name, default = NA) if (name %in% colnames(pd)) pd[[name]] else rep(default, nrow(pd))
cells <- sub("\\.bam$", "", as.character(column("name", rownames(pd))))
out <- data.frame(cell           = cells,
                  ploidy         = as.numeric(column("ploidy")),
                  rpc            = as.numeric(column("rpc")),
                  used_reads     = as.numeric(column("used.reads")),
                  failure_reason = as.character(column("failure_reason")),
                  stringsAsFactors = FALSE)
# A cell that failed to scale carries a meaningless ploidy; blank it so that it is counted as
# unusable downstream instead of being scored as a wrong answer.
failed <- !is.na(out$failure_reason) | !is.finite(out$rpc) | out$rpc <= 0
out$ploidy[failed] <- NA
write.table(out, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("scAbsolute ploidy calls written to:", output_tsv,
    "(", sum(!failed), "of", nrow(out), "cells succeeded )\n")
