#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
from math import log2
import argparse

"""Convert log2FC to FC and gene_id to gene_name. Flag genes signif. up or down vs reference level"""

def get_args():
    """Get arguments from command line"""

    parser = argparse.ArgumentParser(
              description='Annotate the DESeq2 reuslts files')
    parser.add_argument('files', help='file-of-files of DESeq2 output csv files, one per line')
    parser.add_argument('genes', help='gene_id to gene_name conversion table')
    parser.add_argument('-f','--foldchange', help='fold-change cutoff', type=float, default=2.0)
    parser.add_argument('-p', '--padj', help='adjusted p-value cutoff', type=float, default=0.05)
    parser.add_argument('-t', '--tag', type=str, default=str(''),
        help='optional tag to add to output filenames e.g. non-default parameters')
    arguments = parser.parse_args()

    return arguments


def load_conversion(gene_table):
    # read in conversion table
    tbl = {}
    with open(gene_table) as infile:
        for line in infile:
            gid, gname = line.strip().split(' ')
            tbl[gid] = gname
    
    return tbl


def signif_caller(ser, pcut=0.05, fcut=2.0):
    # ser is a pandas series with 2 elements
    padj, log2FC = ser
    fcut_log = log2(fcut)
    if float(padj) < pcut and abs(float(log2FC)) >= fcut_log:
        res = 'YES'
    else:
        res = 'NO'

    return res


def direction_caller(fold_change):
    fc = float(fold_change)
    if fc > 1:
        res = 'UP'
    elif fc < 1:
        res = 'DOWN'
    else:
        res = 'NO_CHANGE'

    return res


def derive_fname(filename):
    """Get the base of a csv filename"""
    fname = filename.strip().split("/")[-1].split('.csv')[0]
    
    return fname


def transform_down(ser):
    """convert the decimal fraction fold change values for the down genes"""
    fc, direction = ser
    f = float(fc)
    if direction == 'UP':
        res = f
    elif direction == 'DOWN':
        res = 1.0/f
    
    return res


def expander(csv_file, con_table, pcut, fcut, tag):
    """ read in a deseq2 results file and expand its information"""
    # load the file
    print(f'Processing {csv_file}')
    colnames = ["gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
    d1 = pd.read_csv(csv_file, names=colnames, header=0)

    # add values
    d1['gene_name'] = d1.gene_id.apply(lambda x: con_table[x])
    d1['fc'] = d1.log2FoldChange.apply(lambda x: 2.0**x)
    d1['significant'] = d1[['padj','log2FoldChange']].apply(signif_caller, axis=1, args=(pcut, fcut))
    d1['change_direction'] = d1.fc.apply(direction_caller)
    d1['fold_change'] = d1[['fc', 'change_direction']].apply(transform_down, axis=1)
    col_order = ['gene_id', 'gene_name', 'baseMean', 'fold_change','log2FoldChange', 'lfcSE',
                 'stat', 'pvalue', 'padj', 'change_direction', 'significant']
    short_order = ['gene_id', 'gene_name', 'fold_change', 'change_direction', 'padj']
    
    d1 = d1[col_order] # reorder columns
    d1 = d1.sort_values(by='padj')
    d2 = d1[d1.significant.isin(['YES'])]
    d2 = d2[short_order]
    d2 = d2.sort_values(by='padj')


    # construct some names
    fname = derive_fname(csv_file)
    # handle presence of a tag nicely
    if len(tag) > 0:
        outname = ''.join([fname, '_', tag, '_full.xlsx'])
        shortname = ''.join([fname, '_', tag, '_DE_only.xlsx'])
        csvname = ''.join([fname, '_', tag, '_full.csv'])
    else:
        outname = ''.join([fname, tag, '_full.xlsx'])
        shortname = ''.join([fname, tag, '_DE_only.xlsx'])
        csvname = ''.join([fname, tag, 'full.csv'])
    
    # write to file
    print(f'Writing {csvname}')
    d1.to_csv(csvname, index=False)
    print(f'Writing {outname}')
    d1.to_excel(outname, index=False, sheet_name=fname, na_rep='ND')
    print(f'Writing {shortname}')
    d2.to_excel(shortname, index=False, sheet_name=fname, na_rep='ND')


def main():
    # process arguments
    args = get_args()
    fof = args.files
    gtable = args.genes
    pcut = args.padj
    fcut = args.foldchange
    tag = args.tag
    
    trans = load_conversion(gtable)
    with open(fof, 'r') as files:
        for file in files:
            f = file.strip()
            expander(f, trans, pcut, fcut, tag)


if __name__ == '__main__':
    main()
