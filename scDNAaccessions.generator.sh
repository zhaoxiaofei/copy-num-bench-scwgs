#!/usr/bin/env bash

function gen_raw() {
csvformat -T SraRunTable.PMID24360273-SRA091188.csv | tr ' ' '~' | awk '{if (NR==1) {h="#";donor="Donor";note="Note"                } else {h="";donor=substr($27, 0, 3);note="NA"         }; print h$1, $4, $23, $33, $(NF-2), $27, donor, $NF,     $5, note}' | sed 's/Oocyte~product:~Second~poloar~body/PB2/g' | sed 's/Oocyte~product:~First~poloar~body/PB1/g' | sed 's/Oocyte~product:~female~pro-nucleus/FPN/g'
csvformat -T SraRunTable.PMC7923680-PRJNA533595.csv | tr ' ' '~' | awk '{if (NR==1) {h="#";donor="Donor";note="Note";ooc="Oocyte_ID"} else {h="";donor=$18;              note="NA";ooc="NA"}; print h$1, $4, $19, $29, $(NF-4), ooc, donor, $(NF-2), $5, note}' | grep -P '345HS1|AvgSpotLen' | sed 's/single~cell/sperm/g' | sed 's/donor~1/345HS1/g' | sort -V -k 2,2 -k 3,3
}

gen_raw | head -n -1 | tr ' ' '\t' # > scDNAaccessions.tsv # generate scDNAaccessions.tsv if the stdout redirection is not commented-out

