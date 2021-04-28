#!/usr/bin/env Rscript

suppressMessages(require(optparse))
# setup arg parser
option_list <- list(
    make_option(c('-d', '--datadir'), type='character', default=NULL,
        help='data directory with htseq-count output files specified in samples.txt',
        metavar='character'),
    make_option(c('-p', '--padj'), type='numeric', default=0.05,
        help='padj threshold', metavar='numeric'),
    make_option(c('-f', '--foldchange'), type='numeric', default=2,
         help='fold-change threshold', metavar='numeric'),
    make_option(c('-n', '--name'), type='character', default='experiment',
         help='name of comparison being plotted', metavar='character'),
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
if (is.null(opt$datadir)){
  print_help(opt_parser)
  stop("The data directory (-d) must be specified", call.=FALSE)
}

# load the rest of the libs
suppressMessages(require(DESeq2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(pheatmap))

# get args
fc_cut <- opt$foldchange
padj_cut <- opt$padj
count_cut <- opt$mincounts
datadir <- opt$datadir
name <- opt$name
samp <- opt$samples
cond <- opt$condition

# MAIN
## process arguments
log2fc_cut <- log2(fc_cut)
samples <- read.table(samp, header=TRUE)
files <- file.path(datadir,samples$run)
names(files) <- samples$library
condition <- read.table(cond, header=TRUE)
treat <- condition$condition

sampleTable <- data.frame(sampleName = samples, fileName = files, condition = treat)
sampleTable$condition <- factor(sampleTable$condition)

# print to console to catch file sorting errors
print("sampleTable")
print(sampleTable)

## create deseq dataset
# import function doesn't seem to accept string for formula
# so, unfortunately, this needs to be hard-coded every time...
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = datadir, design = ~condition) ### CUSTOMIZE
dds <- dds[rowSums(counts(dds)) >= count_cut,] # drop genes where sum counts is < count_cut
dds$condition <- relevel(dds$condition, ref="control") ### CUSTOMIZE

## perform EDA
vsd <- vst(dds, blind=FALSE) # variance stabilizing transform
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(name, "_sample_distances.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf(paste0(name, "_PCA.pdf"))
plotPCA(vsd, intgroup=c("condition")) # CUSTOMIZE
dev.off()

## perform DE according to the specified design
dds <- DESeq(dds)
res <- results(dds, alpha=padj_cut, lfcThreshold=log2fc_cut, altHypothesis='greaterAbs')
write.csv(as.data.frame(res), paste0(name, "_deseq2.csv"))

# write normalized counts to file
write.csv(as.data.frame(counts(dds, normalized=T)), paste0(name, "_deseq2_normCounts.csv"))
