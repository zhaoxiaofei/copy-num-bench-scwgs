#!/usr/bin/env bash

script_dir=$(dirname "$(realpath "$0")")
rootdir=${script_dir}/../

mkdir -p pushd ${rootdir}/refs && pushd ${rootdir}/refs
wget -c https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig
wget -c https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
#wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

tar -vxf chromFa.tar.gz
rm hg19.fa || true
for i in $(seq 1 22) X Y; do cat chr$i.fa >> hg19.fa; done
bwa index hg19.fa
#bowtie2-build hg19.fa hg19.fa
samtools faidx hg19.fa
samtools dict hg19.fa > hg19.fa.dict
cp hg19.fa.dict hg19.dict
popd

ref=${rootdir}/refs/hg19.fa
window_size=200000
bigwig=${rootdir}/refs/wgEncodeCrgMapabilityAlign36mer.bigWig
chrs='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'

time -p mapCounter -w ${window_size} ${bigwig} -c ${chrs} > ${ref}.mp.seg     #parallel=setup.ref/
time -p gcCounter  -w ${window_size} ${ref}    -c ${chrs} > ${ref}.gc.seg     #parallel=setup.ref/
time -p bedtools makewindows -g ${ref}.fai -w 250000      > ${ref}.250000.bed #parallel=setup.ref/
time -p bedtools makewindows -g ${ref}.fai -w 500000      > ${ref}.500000.bed #parallel=setup.ref/

