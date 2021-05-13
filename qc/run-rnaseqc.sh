#!/bin/bash

aln=$1; shift
gtf=$1; shift
outdir=$1

RNASEQC=rnaseqc.v2.3.6.linux # modify as needed

$RNASEQC $gtf $aln $outdir
