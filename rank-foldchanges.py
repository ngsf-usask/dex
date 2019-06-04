#!/usr/bin/env python3
# coding: utf-8

"""for enrichr, weigh the data by fold change and p-value, and scale 0 - 1"""

import pandas as pd
import argparse


def get_args():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
              description='Rank DE genes with function FC/padj')
    parser.add_argument('file', help='DESeq2 output csv file')
    parser.add_argument('-t', '--tag', type=str, default=str(''),
        help='optional tag to add to output filenames e.g. non-default parameters')
    arguments = parser.parse_args()

    return arguments


def derive_fname(filename):
    """Get the base of a csv filename"""
    fname = filename.strip().split("/")[-1].split('.csv')[0]
    
    return fname


def load_data(datfile):
    """load data from a processed deseq2 results csv"""
    cols = ['gene_id','gene_name','baseMean','fold_change','log2FoldChange','lfcSE',
        'stat','pvalue','padj','change_direction','significant']
    dat = pd.read_csv(datfile,header=None, names = cols)
    dat2 = dat[dat.significant == 'YES']

    return dat2


def main():
    args = get_args()
    datfile = args.file
    tag = args.tag
    dat = load_data(datfile)

    minp = dat[dat.padj != 0].padj.min() * 0.01
    dat['padj'].replace(0.0, minp)
    dat['sf'] = dat.apply(lambda x: (x['fold_change'] / x['padj']), axis=1)
    dat.sort_values(by='sf', inplace=True, ascending=False)
    dat['rank'] = dat.sf.rank(method='first', pct=True)

    # handle presence of a tag nicely
    fname = derive_fname(datfile)
    if len(tag) > 0:
        outname = ''.join([fname, '_', tag, '_ranked.csv'])
    else:
        outname = ''.join([fname, tag, '_ranked.csv'])

    dat.to_csv(outname, na_rep=1, index=False, header=False,
        columns=['gene_name', 'rank'])


if __name__ == '__main__':
    main()
