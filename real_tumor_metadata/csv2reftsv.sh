#/stor/zxf/cnv/copy-num-bench-scwgs/real_tumor/SRP002535_Navin100_SraRunTable.csv
in=$1 # ${1/%.sh/}
out=${1}_ref.tsv
cat ${in} | csvcut -c 'Run','AvgSpotLen','BioProject','Library Name','Sample Name','BioSample','SRA Study','Bases','Bytes','LibrarySource','Platform','LibraryLayout' | csvformat -T | awk -F'\t' 'BEGIN{OFS=FS} {if(NR==1) print $0, "Donor"; else print $0, $5}' > ${out}

