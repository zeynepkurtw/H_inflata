import pandas as pd
from Bio import SeqIO

og_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt"
og_count_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
singleton_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups_UnassignedGenes.tsv"

def get_og_stack() -> pd.DataFrame:
    """
    Get og with gene list.
    :input: og.csv
    @return: dataframe with stack genes and their corresponding og
    """

    og_gene_list = pd.read_csv(og_file, header=None, dtype=str, delim_whitespace=True)
    og_gene_list = og_gene_list.apply(lambda x: x.str.replace(":", ""))
    df = og_gene_list.rename(columns={0: "OG"}).set_index("OG")
    df = df.stack().reset_index().drop(columns=["level_1"])
    return df

def get_og_genes(og_category: pd.DataFrame) -> pd.DataFrame:
    """

    @param og_category: Orthgroups with a species combination.
    @return: Orthogroup with its genes for given category.
    """
    df = pd.merge(get_og_stack(), og_category, on="OG")
    return df

def get_og_count() -> pd.DataFrame:
    """
    @return: Orthogroup count file. Order species. Rename orthogroup column.
    """
    df = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df = df[["OG", "carpe", "kbiala", "HIN", "trepo", "spiro", "wb", "muris", "Total"]]
    return df

def get_og_count_sing() -> pd.DataFrame:
    """
    @return: merge og count with singletons
    """
    # OG with at least two genes
    df_count = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df_count = df_count.set_index("OG").sort_values(by="Total", ascending=False)
    df_count.loc[df_count["Total"] > 1, "Type"] = "OG"

    #"Single genes which are excluded from the OG"
    df_sing = pd.read_csv(singleton_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df_sing = df_sing.set_index("OG").fillna(0)
    df_sing["Total"] = 1
    df_sing = df_sing.applymap(lambda x: 1 if isinstance(x, str) == True else x)

    #"Concatanate OG and singleton dataframes"
    df_count_s = pd.concat([df_count, df_sing], axis=0)
    df_count_s.loc[df_count_s["Total"] > 1, "Type"] = "OG"
    df_count_s.loc[df_count_s["Total"] == 1, "Type"] = "singleton"

    return df_count_s.reset_index()

def get_annot(fasta_file, interpro_file):
