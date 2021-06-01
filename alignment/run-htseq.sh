#!/bin/bash

# bam must be position-sorted

label=$1; shift
infile=$1; shift
gtf=$1

htseq-count -f bam \
    -r pos \
    -s reverse \
    -m union \
    --nonunique all \
    $infile \
    $gtf > ${label}.counts.htseq.txt
