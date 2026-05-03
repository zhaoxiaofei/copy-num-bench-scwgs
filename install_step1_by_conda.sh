#!/usr/bin/env bash

if [ -z "${1}" ] ; then
    envname=copy-num-bench-scwgs;
else 
    envname=$1;
fi
envparam="-n ${envname}"

bioconda_packages="hmmcopy bioconductor-hmmcopy bioconductor-ctc bioconductor-dnacopy bioconductor-copynumber bioconductor-wgsmapp bioconductor-scope bioconductor-helloranges"
condaforge_packages="apache-airflow docopt picard pyfaidx r-devtools r-inline r-gplots r-scales r-plyr r-ggplot2 r-gridExtra r-fastcluster r-heatmap3"

conda=mamba
condaforge="" #" -c conda-forge " # conda-forge
bioconda="" #"-c conda-forge -c bioconda"

$conda create -n $envname "python<3.13"  matplotlib numpy 'pandas<2.0' scipy scikit-learn seaborn tasklogger vim bcftools bedtools bowtie2 bwa gatk4 htslib samtools htslib pysam ${condaforge_packages} ${bioconda_packages}

# Install commonly used bioinformatics tools (already installed during $envname creation)
#$conda install --yes $envparam $bioconda bcftools bedtools bowtie2 bwa gatk4 htslib samtools htslib # pysam 

# Install HMMcopy CopyNumber (already installed during $envname creation)
#$conda install --yes $envparam $bioconda bioconductor-hmmcopy hmmcopy bioconductor-copynumber

# Install https://github.com/deepomicslab/SeCNV (already installed during $envname creation)
#$conda install --yes $envparam docopt picard pyfaidx

# Install ginkgo (already installed during $envname creation)
#$conda install --yes $envparam \
#    bioconductor-ctc bioconductor-dnacopy bioconductor-wgsmapp bioconductor-scope bioconductor-helloranges \
#    r-inline r-gplots r-scales r-plyr r-ggplot2 r-gridExtra r-fastcluster r-heatmap3 r-devtools

# Install https://github.com/biosinodx/SCCNV
# echo nothing new to be installed for SCCNV

# Install SCOPE (please install with R if installation error occurs with conda)
# echo nothing new to be installed for SCOPE # $conda install --yes $envparam $bioconda bioconductor-scope 

# Install https://github.com/xikanfeng2/SCYN
$conda install -n $envname --yes tasklogger
$conda run     -n $envname       pip install scyn

# Please uncomment-out if such error occurs on some older versions of pandas
#sed -i 's;all_cnv = all_cnv.append(cnv);all_cnv = pd.concat([all_cnv, cnv]);g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py
#sed -i 's;all_cnv = all_cnv.append(cnv);all_cnv = pd.concat([all_cnv, cnv]);g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py

# Install AneuFinder from Bioconductor
$conda install -y bioconductor-aneufinder

# Install FLCNA dependencies and FLCNA from GitHub
$conda install -y r-base r-devtools bioconductor-genomicranges r-mclust
Rscript -e 'library(devtools); install_github("FeifeiXiao-lab/FLCNA")'

# Install third-party packages to construct the reference genome aux files for FLCNA
$conda install r-essentials \
   bioconductor-bsgenome.hsapiens.ucsc.hg19 \
   bioconductor-genomicranges \
   bioconductor-genomeinfodb \
   bioconductor-iranges \
   bioconductor-annotationhub \
   bioconductor-rtracklayer

# Chisel from https://github.com/raphael-group/chisel?tab=readme-ov-file#automatic
rm chisel || true
git clone https://github.com/raphael-group/chisel && pushd chisel # Clone CHISEL and enters the directory
bash install_full.sh # Run the automatic installation, which installs conda and CHISEL on it
popd

# Alleloscope (encountered running error, as described at https://github.com/seasoncloud/Alleloscope/issues/16)
# conda create -n alleloscope -y r-base r-devtools bioconductor-rtracklayer r-pak r-data.table r-ggplot2 r-remotes r-matrix samtools bcftools r-rlang
# conda run -n alleloscope Rscript -e 'devtools::install_github("seasoncloud/Alleloscope")'

if false; then
# Add alleloscope after chisel when the issue on alleloscope will be solved
for en in ${envname} chisel; do
    conda list       -n ${en} -e | grep -v "^scyn="   > env/${en}.requirements.list_e_no_pypi.txt
    conda env export -n ${en}                         > env/${en}.freeze.env_export.yml
    conda env export -n ${en} --no-builds             > env/${en}.freeze.env_export_no_builds.yml
    conda env export -n ${en} --from-history          > env/${en}.freeze.env_export_from_history.yml
done
fi

# Fix error in scyn and SCOPE due to zero samples remaining after filtering
conda activate $envname && sed -i 's;Gini<=0.12;Gini<=0.21;g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py
conda activate $envname && sed -i 's;perform_qc(Y_raw = Y_raw,;perform_qc(mapq20_thresh = 0.1, Y_raw = Y_raw,;g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py

# The following commented-out code how to install some packages via other package managers
# (such as bioconductor or other conda channels)

#script_dir=$(dirname "$(realpath "$0")")
#pushd "${script_dir}/../software"

#HMM_PATH=/envs/test/lib/R/library/HMMcopy/doc/HMMcopy.R
#conda install -c shahcompbio pypeliner
#git clone https://github.com/KChen-lab/MEDALT.git # 3d4a6d548171ede333310d2ef25c12cdccd11a2b
#Rscript -e 'install.packages("igraph")'
#Rscript -e 'BiocManager::install("HelloRanges")'

#wget https://github.com/shahcompbio/single_cell_pipeline/archive/refs/tags/v0.8.26.tar.gz
#tar -xvf v0.8.26.tar.gz
#rm 1.0.tar.gz || true
## MEDALT performs clustering based on CNV-calling results, instead of calling CNV by itself
#wget https://github.com/KChen-lab/MEDALT/archive/refs/tags/1.0.tar.gz
#tar -xvf 1.0.tar.gz
#git clone https://github.com/KChen-lab/MEDALT.git # 3d4a6d548171ede333310d2ef25c12cdccd11a2b
#Rscript -e 'install.packages("igraph")'
#Rscript -e 'BiocManager::install("HelloRanges")'
#popd

