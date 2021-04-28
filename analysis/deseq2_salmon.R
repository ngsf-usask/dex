#!/usr/bin/env Rscript

suppressMessages(require(optparse))
# setup arg parser
option_list <- list(
    make_option(c('-d', '--datadir'), type='character', default=NULL,
        help='data directory with subdirs named as in samples.txt, each containing quant.sf',
        metavar='character'),
    make_option(c('-p', '--padj'), type='numeric', default=0.05,
        help='padj threshold', metavar='numeric'),
    make_option(c('-f', '--foldchange'), type='numeric', default=2,
         help='fold-change threshold', metavar='numeric'),
    make_option(c('-n', '--name'), type='character', default='experiment',
         help='name of comparison being plotted', metavar='character'),
    make_option(c('-t', '--tx2gene'), type='character', default=NULL,
         help='full path to tx2gene table', metavar='character'),
    make_option(c('-c', '--mincounts'), type='numeric', default=10,
         help='minimum counts per gene per row', metavar='numeric'),
    make_option(c('-s', '--samples'), type='character', default=NULL,
         help='full path to samples.txt', metavar='character'),
    make_option(c('-o', '--condition'), type='character', default=NULL,
         help='full path to file containing condition table', metavar='character')
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check args
if (is.null(opt$datadir)||is.null(opt$tx2gene)){
  print_help(opt_parser)
  stop("The data directory (-d) and tx2gene table (-t) must be specified", call.=FALSE)
}

# load the rest of the libs
suppressMessages(require(DESeq2))
suppressMessages(require(readr))
suppressMessages(require(tximport))
suppressMessages(require(ggplot2))

# get args
fc_cut <- opt$foldchange
padj_cut <- opt$padj
count_cut <- opt$mincounts
datadir <- opt$datadir
tx2geneFile <- opt$tx2gene
tx2gene <- read_csv(tx2geneFile)
name <- opt$name
samp <- opt$samples
cond <- opt$condition

# MAIN
log2fc_cut <- log2(fc_cut)
samples <- read.table(samp, header=TRUE)
files <- file.path(datadir,samples$run,"quant.sf")
names(files) <- samples$library
# the gencode transcript fasta has crazy long header lines, hence ignoreAfterBar
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)
sampleTable <- read.csv(cond, colClasses='factor')
rownames(sampleTable) <- colnames(txi$counts)
# import function doesn't seem to accept string for formula
# so, unfortunately, this needs to be hard-coded every time...
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~clone) #TODO MODIFY FORMULA HERE
dds <- dds[rowSums(counts(dds)) >= count_cut,] # drop genes where sum counts is < count_cut

# perform DE according to the specified design
dds <- DESeq(dds)
# the contrasts argument also needs to be changed PER EXPERIMENT
res <- results(dds, alpha=padj_cut, lfcThreshold=log2fc_cut, altHypothesis='greaterAbs', contrast=c(0,0.5,0.5))
write.csv(as.data.frame(res), paste0(name, "_deseq2.csv"))
