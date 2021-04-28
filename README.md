# dex
Differential expression workflow for RNA-seq data

0. preprocessing/run-fastp.sh

    gently trim reads

1a. mapping/run-salmon.sh

   do quant

1b. alignment/sbatch-star and alignment/run-htseq.sh

   align and quant

2. analysis/deseq2_salmon.R or analysis/deseq2_htseq.R

   do DGE

3. analysis/process-results.py

   annotate genes

4. analysis/interactive-volcano.R

   make interactive and static plots  
   uses cropper.sh and volcano-functions.R

5. analysis/rank-foldchanges.py

   rank FCs for enrichr pathway analysis

6. analysis/IPA_prep.sh

   make input file for IPA

7. analysis/plot-venn.R

   make 2-Venn diagram

