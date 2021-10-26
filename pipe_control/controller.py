#!/usr/bin/env python3

import argparse
import subprocess
import time
import checks

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

list_of_jobs = [] # Keeps all job information in dictionaries


# TODO : build a function to check the arguments

def get_args():
    """
    Purpose:
        Obtain arguments from command line
    Return:
        Object where parameters contain string information from command line.
            .libraries = file name for .txt containing library IDs
            .genomics = file name for .txt containing following information:
                        [Line 1]: Absolute path to raw data directory
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
        subprocess.run(["sbatch", "./batch_pipe_RNAseq.sh",  library.strip(), paths[0]])


    time.sleep(5) # Give nodes enough time to start
    jobIDs = get_jobIDs()
    build_job_list(jobIDs)
    
    # TODO SHOULD I RUN THIS ASYNC?
    # TODO convert this into a function to run ASYNC?
    for library in list_of_jobs:
        started = check_for_start(library)
        print(f"Run has started for {library['lib_ID']}")
        while not(check_for_combine(library)):
            print(f"Checking for when lane combining has finished for {library['lib_ID']}")
            time.sleep(60)
            # TODO is waiting for time the best approach here?


def build_job_list(jobs):
    """
    Purpose:
        Build a library dictionary for every job assigned to plato.
    Pre-cond:
        :param jobs: List of jobIDs provided by get_jobIDs() function
    Return:
        list_of_jobs will be appended to contain all library information
    """
    for jobID in jobs:
        job_info = {}
        job_info["jobID"] = jobID
        job_info["slurm"] = f"slurm-{jobID}.out"
        job_info["lib_ID"] = ""
        # Booleans to track when steps are complete
        job_info["start"] = False
        job_info["done_combining"] = False
        job_info["read1_counts"] = 0
        job_info["read2_counts"] = 0
        job_info["done_fastp"] = False
        job_info["done_STAR"] = False
        job_info["done_sbatch_job1"] = False
        list_of_jobs.append(job_info)

def get_jobIDs():
    """
    Purpose:
        Obtain a list of all active jobIDs from "NGSF_RNA" runs on sbatch nodes
    Return:
        List of jobIDs active on the compute cluster
    """
    queue = subprocess.run(["squeue", "-n", "NGSF_RNA"], capture_output=True)
    jobIDs = []
    for line in queue.stdout.splitlines():
        steps = (str(line).split())
        if steps[1].isdigit():
            jobIDs.append(str(steps[1]))
    return jobIDs


def cancel_all_runs():
    """
    Purpose:
        If the RNA-seq fails, this will cancel every batch run.
    Return:
        Cancels every batch run.
        Kills python script.
        - Untested
    """
    for job in list_of_jobs:
        subprocess.run(["scancel", job["jobID"]])
    exit()


def check_for_start(library):
    """
    Purpose:
        Checks to ensure node started correctly
    Returns:
        True if node started correctly.
        Cancels script otherwise.
    """
    if checks.start(library["slurm"]):
        library["start"] = True
        library["lib_ID"] = checks.library(library["slurm"])
        return True
    else:
        print(f"Library failed to start properly.")
        cancel_all_runs()

def check_for_combine(library):
    """
    Purpose:
        Checks to four lanes of data were combined into one file. 
    Returns:
        True if data has been combined
        False otherwise.
    """
    # TODO how do I check to see if combine failed?
    if checks.combined(library["slurm"]):
        library["read1_counts"] = checks.read1(library["slurm"])
        library["read2_counts"] = checks.read2(library["slurm"])
        return True
    else:
        return False

def check_for_completion(library):
    started = checks.start(library["slurm"])
    if started:
        library["start"] = True
    else:
        pass
        # print(f"{library[")

def main():
    args = get_args()
    call_batch_runs(args.libraries, args.genomics)

if __name__ == "__main__":
    main()