#!/usr/bin/env R

suppressMessages(require(readr))
suppressMessages(require(rbokeh))
suppressMessages(require(htmlwidgets))
suppressMessages(require(webshot))

load_data <- function(csv){
    col_types <- cols(
        gene_id = col_character(),
        gene_name = col_character(),
        baseMean = col_double(),
        fold_change = col_double(),
        log2FoldChange = col_double(),
        lfcSE = col_double(),
        stat = col_double(),
        pvalue = col_double(),
        padj = col_double(),
        change_direction = col_character(),
        significant = col_character()
    )

    dat <- read_csv(datf, col_types=col_types)

    # set padj == 0 to 1/10 the lowest padj for plotting
    minp <- min(dat$padj[which(dat$padj>0)]) * 0.1
    dat$padj <- replace(dat$padj, dat$padj == 0, minp)
    dat$padj <- replace(dat$padj, is.na(dat$padj), 1)

    dat['-log10padj'] <- -1*log(dat$padj)
    dat['Significant'] <- apply(dat, 1, signif_caller, pcut=pcut, fcut=fcut)
    dat['color'] <- apply(dat, 1, color_selector)
    datS <- dat[which(!(dat$Significant == 'Not Significant')),]

    res <- list(dat, datS)
    names(res) <- c('full', 'signif')

    return(res)
} 
        

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


shooter <- function(webfile_name, tries=10){
    # webfile_name is partial filename of html rbokeh plot to capture
    message('Trying to capture web plot')
    prober <- function(webfile_name, trial){
        if(trial==1){
            message(paste0('Try ', trial))
        }
        else{
            message(paste0('\nTry ', trial))
        }

        out <- tryCatch(
            {webshot(paste0(name, '_volcano_full.html'), paste0(name, '_volcano_fullX.png'))},
            error=function(e){
                message('Failed to capture volcano plot. Error message:')
                message(e)
                return(NA)
            }
            )
            return(out)
        }
    for(trial in 1:tries){
        result <- prober(webfile_name, trial)
        if(!(is.na(result))){
            message('Web plot successfully captured!')
            break
        }
    }
}
