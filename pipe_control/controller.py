#!/usr/bin/env python3

import argparse
import os
import subprocess

"""
This script will control the data pipeline for RNA seq data.
Steps.
1. Process input .txt files containing library information, and pathways for indexed genomes
2. Run batch-RNA-seq.sh on a plato node
    -> Data from 4 lanes will be combined into a single file.
    -> Run fastp to QC to reads
        # Report is exported from the node here
    -> Run STAR using provided genome
        # Output all STAR documents
3. The processed results will be analyzed using process_results.py on a node
4. Volcano plots, and any other results will be exported out of the node.
"""

def get_args():
    """
    Purpose:
        Obtain arguments from command line
    Return:
        Object where parameters contain string information from command line.
            .libraries = file name containing library IDs
            .genomics = file name containing path to genome
    """
    parser = argparse.ArgumentParser(
        description="Process raw nextSeq RNA-seq results")
    parser.add_argument("libraries", help="txt file where each line is a library ID")
    parser.add_argument("genomics", help="txt file where first line is absolute path to raw data")
    arguments = parser.parse_args()
    return arguments

def call_batch_runs(lib_file, genome_file):
    gen_file = open(genome_file, "r")
    paths =[]
    for path in gen_file:
        paths.append(path.strip())

    in_file = open(lib_file, "r")
    for library in in_file:
        os.system(f"sbatch ./batch_pipe_RNAseq.sh {library.strip()} {paths[0]}")
        # TODO keep working on batch pipe



def main():
    args = get_args()
    call_batch_runs(args.libraries, args.genomics)
    print(args.libraries)


if __name__ == "__main__":
    main()