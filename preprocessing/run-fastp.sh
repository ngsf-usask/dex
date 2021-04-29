#!/bin/bash

# run fastp to trim adapters and low-quality sequences from reads
# usage run-fastp.sh file-of-files.txt
# file-of-files has R1 and R2 on single line separated by a space

THREADS=6
infile=$1

while read line; do
    file1=$(echo $line | awk '{print $1}')
    fname1=$(echo $file1 | awk -F\/ '{print $NF}')
    bname1=$(echo $fname1 | sed 's/.fastq.gz//')
    file2=$(echo $line | awk '{print $2}')
    fname2=$(echo $file2 | awk -F\/ '{print $NF}')
    bname2=$(echo $fname2 | sed 's/.fastq.gz//')
    corename=$(echo $bname1 | sed 's/_R1//')

    $fastp $line --in1 $file1 --in2 $file2 \
        --out1 ${bname1}_trimmed.fastq.gz --out2 ${bname2}_trimmed.fastq.gz \
        --unpaired1 ${bname1}_trimmed_unpaired.fastq.gz \
        --unpaired2 ${bname2}_trimmed_unpaired.fastq.gz \
        -V \
        -l 30 \
        -p \
        -w $THREADS \
        -j ${corename}.trim-report.json \
        -h ${corename}.trim-report.html \
        -e 10 \
        -q 10 \
        -M 10 \
        -r \
        -W 6 \
        -g
done < $infile
