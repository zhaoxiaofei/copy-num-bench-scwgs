#!/usr/bin/env bash
set -euo pipefail

# scAbsolute -- only needed by the opt-in `--ploidy-tools scabsolute` benchmark.
# Run explicitly from this directory:
#   bash -evx install_scabsolute.sh
#
# The authoritative environment is the lab's own pipeline:
#   https://github.com/markowetzlab/scDNAseq-workflow
#
# This creates a self-contained conda environment containing both the R and
# Python dependencies needed by scAbsolute.

# ----------------------------------------------------------------------
# Create environment
# ----------------------------------------------------------------------

conda create --yes --name scabsolute \
    --override-channels \
    -c conda-forge \
    -c bioconda \
    python \
    r-base \
    r-reticulate \
    r-future.apply \
    r-tidyverse \
    r-devtools \
    r-digest \
    r-robustbase \
    r-matrixstats \
    r-mass \
    bioconductor-qdnaseq \
    bioconductor-biobase \
    bioconductor-genomicranges \
    bioconductor-rsamtools \
    bioconductor-qdnaseq.hg19 \
    numpy \
    pandas \
    scipy \
    "tensorflow=2.6" \
    "tensorflow-probability<0.15" \
    c-compiler \
    make

# ----------------------------------------------------------------------
# Clone scAbsolute
# ----------------------------------------------------------------------

if [[ ! -d scAbsolute ]]; then
    git clone https://github.com/markowetzlab/scAbsolute.git
fi

# ----------------------------------------------------------------------
# Build scAbsolute's bundled changepoint C code
# ----------------------------------------------------------------------
#
# Compile with the R/compiler toolchain from the scabsolute environment,
# rather than with R from the caller's currently active environment.

pushd scAbsolute/data/changepoint

conda run -n scabsolute \
    R CMD SHLIB \
    cost_general_functions.c \
    PELT_one_func_minseglen.c \
    -o PELT.so

popd

# ----------------------------------------------------------------------
# Verify Python dependencies
# ----------------------------------------------------------------------
#
# scAbsolute's scaling step runs through reticulate and imports numpy,
# pandas, tensorflow and tensorflow_probability.

conda run -n scabsolute python - <<'PY'
import numpy
import pandas
import scipy
import tensorflow as tf
import tensorflow_probability as tfp

print("Python deps OK")
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("scipy:", scipy.__version__)
print("tensorflow:", tf.__version__)
print("tensorflow_probability:", tfp.__version__)
PY

# ----------------------------------------------------------------------
# Verify QDNAseq annotations
# ----------------------------------------------------------------------
#
# QDNAseq only ships bin annotations for a fixed set of bin sizes.
# 500 kb (the working point of the scAbsolute paper, and the default of
# simplerun_scabsolute.R) is one of them.

conda run -n scabsolute \
    Rscript -e 'QDNAseq::getBinAnnotations(binSize=500, genome="hg19")'

# ----------------------------------------------------------------------
# Optional: point ploidy_tools.py at this clone
# ----------------------------------------------------------------------
#
# If scAbsolute was not cloned next to the calling script:
#
#   export scAbsoluteRoot=/path/to/scAbsolute
