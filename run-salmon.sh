#!/bin/bash

r1=$1;shift
r2=$1; shift
prefix=$1; shift
index=$1
#index=/mnt/NGSF001/analysis/indices/human/gencode-30/salmon/gencode.v30.transcripts.index

/home/ahammond/src/salmon-latest_linux_x86_64/bin/salmon quant -i $index -l A -1 $r1 -2 $r2 \
	--validateMappings \
	-o ${prefix}_salmon_quant \
	-p 6 \
	--rangeFactorizationBins 4 \
	--numBootstraps 10 \
    --seqBias \
    --gcBias \
    --biasSpeedSamp=5
