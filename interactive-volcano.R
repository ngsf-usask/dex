#!/usr/bin/env Rscript

# produce an interactive volcano plot of the DE genes, 
# and a static png of all genes

# arg parsing func
suppressMessages(require(optparse))
getArgs <- function(){
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
    return(opt)
}

# check args
opt <- getArgs()
if (is.null(opt$csv)){
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# load the other functions
filepath <- function(){
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename <- dirname(script.name)

    return(script.basename)
}

script_dir <- filepath()
source(paste0(script_dir, '/volcano-functions.R'))

# collect args
datf <- opt$csv
pcut <- opt$padj
fcut <- opt$foldchange
name <- opt$name

# read the data and add label DE genes
res <- load_data(datf)
full <- res$full
DE <- res$signif

#print(table(full$Significant))
print(table(full$Significant))

# the full dataset is too large for most machines to handle interactively
# plot just the significant genes instead, datS
suppressWarnings(p <- figure(title=name) %>%
    ly_points(log2FoldChange, '-log10padj', data = DE,
    color=color,
    hover = list(gene_name, fold_change, change_direction, padj)))
saveWidget(p, paste0(name, '_volcano.html'))

# plot the full set as png
suppressWarnings(q <- figure(title=name) %>%
    ly_points(log2FoldChange, '-log10padj', data = full,
    color=color
    ))
saveWidget(q, paste0(name, '_volcano_full.html'))
shooter(name) ### TODO finish this fn # TODO all lines below depend on fixing this one

# crop the png with imagemagick through shell, for reasons
cmd <- paste0(script_dir, '/cropper.sh ', paste0(name, '_volcano_fullX.png'))
system(cmd)

# delete the uncropped png
rmcmd <- paste0('rm ', paste0(name, '_volcano_fullX.png'))
system(rmcmd)
