#!/usr/bin/env Rscript

# produce an interactive volcano plot of the DE genes, 
# and a static png of all genes

suppressMessages(require(optparse))
# setup arg parser
option_list <- list(
    make_option(c('-c', '--csv'), type='character', default=NULL,
        help='processed deseq2 results csv', metavar='character'),
    make_option(c('-p', '--padj'), type='numeric', default=0.05,
        help='padj threshold', metavar='numeric'),
    make_option(c('-f', '--foldchange'), type='numeric', default=2,
         help='fold-change threshold', metavar='numeric'),
    make_option(c('-n', '--name'), type='character', default='deseq2',
         help='name of comparison being plotted', metavar='character')
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check args
if (is.null(opt$csv)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# laod the rest of the packages
suppressMessages(require(rbokeh))
suppressMessages(require(readr))
suppressMessages(require(htmlwidgets))
suppressMessages(require(webshot))


signif_caller <- function(tibble_row, pcut=0.05, fcut=2.0){
    # tibble row includes named elements padj and log2FoldChange
    padj = as.numeric(tibble_row['padj'])
    log2FC = as.numeric(tibble_row['log2FoldChange'])
    fcut_log = log2(fcut)
    if(is.na(padj)){
        res = 'Not Significant'
    }
    else if((padj < pcut) & (log2FC >= fcut_log)){
        res = 'Up'
    }
    else if((padj < pcut) & (log2FC <= -1*fcut_log)){
        res = 'Down'
    }
    else{
        res = 'Not Significant'
    }
    return(res)
}


color_selector <- function(tibble_row){
    # tibble row includes named elements padj and log2FoldChange
    colors <- as.list(c('orange', 'blue', 'gray'))
    sig_types <- c('Up', 'Down', 'Not Significant')
    names(colors) <- sig_types
    
    sig = tibble_row['Significant']
    res = as.character(colors[sig])

    return(res)
}


# collect args
datf <- opt$csv
pcut <- opt$padj
fcut <- opt$foldchange
name <- opt$name

# read the data and add label DE genes
dat <- read_csv(datf)
dat['-log10padj'] <- -1*log(dat$padj)
dat['Significant'] <- apply(dat, 1, signif_caller, pcut=pcut, fcut=fcut)
dat['color'] <- apply(dat, 1, color_selector)
datS <- dat[which(!(dat$Significant == 'Not Significant')),]

# the full dataset is too large for most machines to handle interactively
# plot just the significant genes instead, datS
p <- figure(title=name) %>%
  ly_points(log2FoldChange, '-log10padj', data = datS,
    color=color,
    hover = list(gene_name, fold_change, change_direction, padj))
saveWidget(p, paste0(name, '_volcano.html'))

# plot the full set as png
q <- figure(title=name) %>%
  ly_points(log2FoldChange, '-log10padj', data = dat,
    color=color,
    hover = list(gene_name, fold_change, change_direction, padj))
saveWidget(q, paste0(name, '_volcano_full.html'))
webshot(paste0(name, '_volcano_full.html'), paste0(name, '_volcano_fullX.png'))
# crop the png with imagemagick through shell, for reasons
cmd <- paste0('./cropper.sh ', paste0(name, '_volcano_fullX.png'))
system(cmd)
