#!/usr/bin/env Rscript

# Map gene expression changes to a KEGG pathway
# data is 2-column csv with ensembl gene ID (NO revision number) and 
#  fold change with direction indicated
# e.g.
# ENSG00000139278,-12.565519453055266
# ENSG00000138829,3.5477686669692727

# arg parsing func
suppressMessages(require(optparse))
getArgs <- function(){
    # setup arg parser
    option_list <- list(
        make_option(c('-c', '--csv'), type='character', default=NULL,
            help='DE genes, gene_id and fold change with direction', metavar='character'),
        make_option(c('-p', '--pathway'), type='character', default=NULL,
            help='KEGG pathway, e.g. hsa04151', metavar='character'),
        make_option(c('-s', '--species'), type='character', default='hsa',
            help='KEGG species, e.g. hsa', metavar='character'),
        make_option(c('-n', '--name'), type='character', default='kegg',
             help='name of comparison being plotted', metavar='character')
    )

    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

# check args
opt <- getArgs()
if (is.null(opt$csv) || is.null(opt$pathway)){
  print_help(opt_parser)
  stop("Input file and KEGG pathway must be supplied", call.=FALSE)
}

suppressMessages(require(pathview))
suppressMessages(require(readr))

datF <- opt$csv
pway <- opt$pathway
species <- opt$species
tag <- opt$name

datDF <- read_csv(datF, col_names=F)
dat <- datDF$X2
names(dat) <- datDF$X1

pv.out <- pathview(gene.data=dat, pathway.id=pway, species=species,
    gene.idtype=gene.idtype.list[3], kegg.native=T, out.suffix=tag, same.layer=T)
