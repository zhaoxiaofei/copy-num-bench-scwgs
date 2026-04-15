#!/usr/bin/env Rscript
# revised from https://sorryios.ai/chat/e62de0ac-e570-4d62-9672-1ba233c8db90
# simplerun_aneufinder.R
# Usage: Rscript simplerun_aneufinder.R <input_bam_dir> <output_dir> <output_bed>
# Runs AneuFinder on BAM files in input_bam_dir, then decompresses the CNV BED.gz to output_bed

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript simplerun_aneufinder.R <input_bam_dir> <output_dir> <output_bed>")
}

input_dir  <- args[1]
output_dir <- args[2]
output_bed <- args[3]

library(AneuFinder)
# pdf(file = "/dev/null")  # suppress plot output
pdf(file = paste0(output_bed, "_hotspot_plots.pdf"), onefile = TRUE, width = 10, height = 8)
r <- Aneufinder(inputfolder = input_dir, outputfolder = output_dir)
dev.off()

# Decompress the CNV BED.gz output to a flat file for downstream parsing
bed_gz <- Sys.glob(file.path(output_dir, "BROWSERFILES", "method-edivisive", "*_CNV.bed.gz"))
if (length(bed_gz) == 0) {
    stop("AneuFinder did not produce the expected CNV BED.gz output!")
}
# Use the first (and typically only) matching file
con_in  <- gzfile(bed_gz[1], "rt")
lines   <- readLines(con_in)
close(con_in)
writeLines(lines, output_bed)
cat("AneuFinder output written to:", output_bed, "\n")

