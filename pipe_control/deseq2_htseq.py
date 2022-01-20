import pandas as pd
import numpy as np
import os

from py_deseq import py_DESeq2


def create_count_matrix(libraries: list, outdir: str) -> pd.DataFrame:
    """
    Purpose:
        Combines htseq results from all provided libraries in a given directory. 
    Pre:
        :param libraries: List of all libraries to analyze
        :param outdir: Relative pathway to project pathway.
    Return:
        Pandas dataframe containing the gene names and the htseq counts for each library.
    """
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

def deseq2(count_matrix, design_matrix):
    dds = py_DESeq2(count_matrix=count_matrix[count_matrix["R2100207"] != 0], # Remove any that contain 0s
                    design_matrix=design_matrix,
                    design_formula="~ treatment",
                    gene_column="Gene")

    dds.run_deseq()
    dds.get_deseq_result(contrast=["treatment", "control", "treated"])
    res = dds.deseq_result
    print(res.head())

def main():
    count_matrix = create_count_matrix(["R2100207", "R2100208", "R2100209", "R2100210"], ".")
    design_matrix = create_design_matrix([["sample", "treatment", "other"],
                    ["R2100207", "control", "other1"],
                    ["R2100208", "control", "other1"],
                    ["R2100209", "treated", "other2"],
                    ["R2100210", "treated", "other2"]])
    deseq2(count_matrix, design_matrix)
    # print(len(count_matrix.index))
    # # df = count_matrix.loc[((count_matrix!=0).any(axis=1))]
    # # df = count_matrix[count_matrix["R2100207"] != 0]
    # print(len(df.index))


if __name__ == "__main__":
    print("Testing?")
    main()
