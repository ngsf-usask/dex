import pandas as pd
import numpy as np
import os

# from py_deseq import py_DESeq2


def combine_results(libraries: list, outdir: str) -> pd.DataFrame:
    # Create the first dataframe
    # TODO: SHOULD I REMOVE ANY 0s?
    df = pd.read_table(f"{outdir}/alignment/{libraries[0]}/{libraries[0]}.counts.htseq.txt", 
                        names=["Gene", libraries[0]])

    # merge any additional libraries
    if len(libraries) > 1:
        for n in range(1, len(libraries)):
            df_to_add = pd.read_table(f"{outdir}/alignment/{libraries[n]}/{libraries[n]}.counts.htseq.txt",   names=["Gene", libraries[n]])
            df = df.merge(df_to_add, how="outer", on="Gene")

    return df

def create_design_matrix(list_matrix: list) -> pd.DataFrame:
    df = pd.DataFrame(data=list_matrix[1:], columns=list_matrix[0])
    df = df.set_index("sample")
    return df

def main():
    # results = combine_results(["R2100207", "R2100208"], ".")
    design_matrix = create_design_matrix([["sample", "treatment", "other"],
                    ["R2100207", "control", "other1"],
                    ["R2100208", "treated", "other2"]])

if __name__ == "__main__":
    print("Testing?")
    main()
