#!/usr/bin/env Rscript

require(DESeq2)
require(readr)
require(tximport)
require(ggplot2)

# GLOBALS
fc_cut <- 2
padj_cut <- 0.05
count_cut <- 10
datadir <- "/mnt/NGSF001/experiments/anderson/demo_run/expression/gencode-30/"
tx2gene <- read_csv('/mnt/NGSF001/analysis/references/human/gencode-30/gencode-30-tx2gene.csv')

# MAIN
# general setup for 2 conditions, 2 samples per condition
log2fc_cut <- log2(fc_cut)
samples <- read.table(paste0(datadir, "stats/samples.txt"), header=TRUE)
samples <- samples[c(1:4),]
files <- file.path(datadir,samples$run,"quant.sf")
files <- files[c(1:4)]
names(files) <- samples$library
# the gencode transcript fasta has crazy long header lines
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)
sampleTable <- data.frame(condition=factor(rep(c('A', 'B'), each=2)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- dds[rowSums(counts(dds)) >= count_cut,] # drop genes where sum counts is < count_cut

# perform DE according to the specified design
dds <- DESeq(dds)
res <- results(dds, alpha=padj_cut, lfcThreshold=log2fc_cut, altHypothesis='greaterAbs')
write.csv(as.data.frame(res), "deseq2.csv")
