# scAbsolute -- only needed by the opt-in `--ploidy-tools scabsolute` benchmark, so this script
# is kept out of install_soft3to4.sh and has to be run explicitly, from this directory:
#   bash -evx install_scabsolute.sh
# The authoritative environment is the one of the lab's own pipeline, and this is a minimal
# stand-in for it: https://github.com/markowetzlab/scDNAseq-workflow

git clone https://github.com/markowetzlab/scAbsolute.git

# scAbsolute segments with its own copy of the changepoint C code, which ships as source.
pushd scAbsolute/data/changepoint
R CMD SHLIB cost_general_functions.c PELT_one_func_minseglen.c -o PELT.so
popd

# The scaling step runs in python through reticulate (tensorflow + tensorflow-probability),
# and QDNAseq supplies the bin annotations, so R and python have to share one environment.
conda create --yes --name scabsolute -c conda-forge -c bioconda \
	python=3.9 r-base=4.2 r-reticulate r-future.apply r-tidyverse r-devtools r-digest \
	r-robustbase r-matrixstats r-mass bioconductor-qdnaseq bioconductor-biobase \
	bioconductor-genomicranges bioconductor-rsamtools bioconductor-qdnaseq.hg19
conda run -n scabsolute pip install "tensorflow>=2.8,<2.16" "tensorflow-probability<0.24"

# QDNAseq only ships bin annotations for a fixed set of bin sizes; 500 kb (the working point of
# the scAbsolute paper, and the default of simplerun_scabsolute.R) is one of them.
conda run -n scabsolute Rscript -e 'QDNAseq::getBinAnnotations(binSize=500, genome="hg19")'

# Point ploidy_tools.py at the clone if it was not put next to this script:
#   export scAbsoluteRoot=/path/to/scAbsolute
