# dex
Differential expression workflow for RNA-seq data

1. run-salmon.sh

   do quant

2. deseq2_base.R

   do DE

3. process-results.py

   annotate genes

4. interactive-volcano.R

   make interactive and static plots
   uses cropper.sh and volcano-functions.R

5. rank-foldchanges.py

   rank FCs for enrichr pathway analysis

6. IPA_prep.sh

   make input file for IPA

7. plot-venn.R

   make 2-Venn diagram

