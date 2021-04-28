#!/usr/bin/env Rscript

# arg parsing func
suppressMessages(require(optparse))
getArgs <- function(){
    # setup arg parser
    option_list <- list(
        make_option(c('-1', '--file1'), type='character', default=NULL,
            help='processed deseq2 results csv', metavar='character'),
        make_option(c('-2', '--file2'), type='character', default=NULL,
            help='processed deseq2 results csv', metavar='character'),
        make_option(c('-n1', '--name1'), type='character', default=NULL,
            help='name of set 1', metavar='character'),
        make_option(c('-n2', '--name2'), type='character', default=NULL,
            help='name of set 2', metavar='character'),
        make_option(c('-n', '--name'), type='character', default='deseq2',
             help='name of comparison being plotted', metavar='character')
    )

    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

# check args
opt <- getArgs()
if (is.null(opt$file1) || is.null(opt$file2)){
  print_help(opt_parser)
  stop("Two input files are required", call.=FALSE)
}

suppressMessages(require(VennDiagram))
suppressMessages(require(readr))

dat1F <- opt$file1
dat2F <- opt$file2
suppressMessages(dat1 <- read_csv(dat1F))
suppressMessages(dat2 <- read_csv(dat2F))

dat1D <- dat1[dat1$significant == 'YES',]$gene_id
dat2D <- dat2[dat2$significant == 'YES',]$gene_id

d12 <- length(intersect(dat1D, dat2D))
d1 <- length(setdiff(dat1D, dat2D))
d2 <- length(setdiff(dat2D, dat1D))

png(paste0(name, '.png'), width=20, height=20, units='cm', res=300)
grid.newpage()
draw.pairwise.venn(d1, d2, d12, category = c(opt$name1, opt$name2), lty = rep("blank", 
    2), fill = c("#008080", "#E34234"), alpha = rep(0.5, 2), cat.pos = c(-8, 8), cat.dist = rep(0.025, 2),
     scaled = TRUE, cex=2, cat.cex=2)
suppressMessages(dev.off())
