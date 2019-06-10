#!/usr/bin/env python3

import pandas as pd
import argparse

""" select DE genes and output gene_id and FC with direction indicated """

def get_args():
    """Get arguments from command line"""

    parser = argparse.ArgumentParser(
              description='Construct input file for IPA')
    parser.add_argument('csv', help='annotated DESeq2 output csv file')
    arguments = parser.parse_args()

    return arguments


def call_dir(x):
    """call direction of a fold change"""
    # x is a row slice
    cd, fc = x
    if cd == 'UP':
        dirFC = fc
    elif cd == 'DOWN':
        dirFC = -1.0 * fc

    return dirFC


def derive_fname(filename):
    """Get the base of a csv filename"""
    fname = filename.strip().split("/")[-1].split('.csv')[0]
    
    return fname


def main():
    # process arguments
    args = get_args()
    infile = args.csv

    colnames = ['gene_id', 'gene_name', 'baseMean', 'fold_change','log2FoldChange', 'lfcSE',
                'stat', 'pvalue', 'padj', 'change_direction', 'significant']
    dat = pd.read_csv(infile, names=colnames, header=0)

    sig = dat[dat.significant.isin(['YES'])]
    sig['dirFC'] = sig[['change_direction', 'fold_change']].apply(call_dir, axis=1)
    sig['gene_id2'] = sig.gene_id.apply(lambda x: x.split('.')[0])

    final = sig[['gene_id2', 'dirFC']]
    
    fname = derive_fname(infile)
    outname = f'{fname}_IPA.tsv'
    final.to_csv(outname, index=False, header=False, sep='\t')


if __name__ == '__main__':
    main()
