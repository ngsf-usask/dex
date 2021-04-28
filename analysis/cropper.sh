#!/bin/bash

# trim the static volcano plot

image=$1
outname=$(echo $image | sed 's/X.png/.png/')

magick $image -crop 470x504+40+40 $outname
