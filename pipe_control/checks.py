"""
Script to control the check steps involved in RNA-seq pipeline.

Each function here takes the slurm.out file and searches for a specific term.
"""

import subprocess

def start(slurm_file):
    """
    Purpose:
        Determines if the run correctly started.
    Return:
        True if slurm output contains a START statement.
    """
    stdout = subprocess.run(["grep", "^#NGSF-START", slurm_file], capture_output=True).stdout.decode("utf-8")
    return True if "True" in stdout else False


def library(slurm_file):
    """
    Purpose:
        Determines the library ID of a given batch run.
    Return:
        Library ID as a string.
    """
    stdout = subprocess.run(["grep", "^#NGSF-LIB", slurm_file], capture_output=True).stdout.splitlines()[0].decode("utf-8").split(" ")[-1]
    return stdout

def combined(slurm_file):
    """
    Purpose:
        Determines if all the four-lane reads have finished combining.
    Return:
        True if slurm output contains a COMBINED statement.
    """
    stdout = subprocess.run(["grep", "^#NGSF-COMBINED", slurm_file], capture_output=True).stdout.decode("utf-8")
    print(stdout)
    return True if "True" in stdout else False

def read1(slurm_file):
    stdout = subprocess.run(["grep", "^#NGSF-R1", slurm_file], capture_output=True).stdout.splitlines()[0].decode("utf-8").split(" ")[-2]
    return stdout

def read2(slurm_file):
    stdout = subprocess.run(["grep", "^#NGSF-R2", slurm_file], capture_output=True).stdout.splitlines()[0].decode("utf-8").split(" ")[-2]
    return stdout