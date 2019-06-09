#!/bin/bash

# select DE genes and output gene_id and FC with direction indicated

# e.g. deseq2_ABvsAA_fc1.4_full.csv 
infile=$1
nam=$(echo $infile | awk -F\/ '{print $NF}' | sed 's/_full.csv//')

awk -F\, '$NF=="YES" && $(NF-1)=="UP"{print $1,$4};$NF=="YES" && $(NF-1)=="DOWN"{print $1,"-"$4}' $infile | sed 's/\./ /' | awk '{print $1","$3}' > ${nam}_dir.csv
