import pandas as pd
from Bio import SeqIO
import glob

og_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt"
og_count_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
singleton_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups_UnassignedGenes.tsv"

fasta_ann_cat = "data/fasta_ann/ann_f_cat.csv"
interpro_ann_cat = "data/interpro_ann/ann_ipr_cat.csv"
eggnog_ann_cat = "data/eggnog_ann/ann_egg_cat.csv"

trepo_file = "data/LGT/trepo_lgt.csv"


def get_og_stack(og_file: str) -> pd.DataFrame:
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


def get_og_genes(og_file, og_category: pd.DataFrame) -> pd.DataFrame:
    """
    @param og_category: Orthgroups with a species combination.
    @return: Orthogroup with its genes for given category.
    """
    df = pd.merge(get_og_stack(og_file), og_category, on="OG")
    return df


def get_og_genes_(og_category: pd.DataFrame) -> pd.DataFrame:
    """
    @param og_category: Orthgroups with a species combination.
    @return: Orthogroup with its genes for given category.
    """
    df = pd.merge(get_og_stack(), og_category, on="OG")
    return df


def get_og_count(og_count_file) -> pd.DataFrame:
    """
    @return: Orthogroup count file. Order species. Rename orthogroup column.
    """
    df = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df = df[["OG", "carpe", "kbiala", "HIN", "trepo", "spiro", "wb", "muris", "Total"]]
    return df


def get_og_count_sing(og_count_file) -> pd.DataFrame:
    """
    @return: merge og count with singletons
    """
    # OG with at least two genes
    df_count = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df_count = df_count[["OG", "carpe", "kbiala", "HIN", "trepo", "spiro", "wb", "muris", "Total"]]
    df_count = df_count.set_index("OG").sort_values(by="Total", ascending=False)
    df_count.loc[df_count["Total"] > 1, "Type"] = "OG"

    # "Single genes which are excluded from the OG"
    df_sing = pd.read_csv(singleton_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
    df_sing = df_sing.set_index("OG").fillna(0)
    df_sing["Total"] = 1
    df_sing = df_sing.applymap(lambda x: 1 if isinstance(x, str) == True else x)

    # "Concatanate OG and singleton dataframes"
    df_count_s = pd.concat([df_count, df_sing], axis=0)
    df_count_s.loc[df_count_s["Total"] > 1, "Type"] = "OG"
    df_count_s.loc[df_count_s["Total"] == 1, "Type"] = "singleton"

    return df_count_s.reset_index()


def get_ann(df):
    """
    Get gene annottaions from fasta and interproscan
    @param df: Genes in the og category to be annotated
    @return: og, genes, ann
    """
    ann_f_cat = pd.read_csv(fasta_ann_cat, sep="\t", header=None).rename(columns={1: "ann_f"})
    ann_ipr_cat = pd.read_csv(interpro_ann_cat, sep="\t", header=None).rename(columns={1: "ipr", 2: "ann_inter"})
    ann_egg_cat = pd.read_csv(eggnog_ann_cat, sep="\t", header=None).rename(
        columns={1: "COG_cat", 2: "KEGG_KOs", 3: "ann_egg"})

    df = pd.merge(df, ann_f_cat, on=0)
    df = pd.merge(df, ann_ipr_cat, on=0)
    df = pd.merge(df, ann_egg_cat, on=0)

    df = df[["OG", "id", "ann_f", "ann_inter", "ann_egg", "ipr",
             "COG", "KEGG_KOs", "carpe", "kbiala", "HIN",
             "trepo", "spiro", "wb", "muris", "Total"]]

    return df


def highlight_rows(df):
    color = (df.LGT == 'trepo').map({True: 'background-color: yellow', False: ''})
    return df.style.apply(lambda s: color)


def mark_trepo_lgt(df):
    trepo = pd.read_csv(trepo_file, sep="\t", header=None)[[0, 1]].rename(columns={0: "OG", 1: 0})
    trepo["LGT"] = "trepo"
    df = pd.merge(df, trepo, how="left").sort_values(by=["LGT", "OG", "Total"], ascending=[False, False, False])
    # highlight_rows(df)
    return df


def merge_trepo_lgt(trepo_file):
    trepo = pd.read_csv(trepo_file, sep="\t", header=None)[[0, 1]].rename(columns={0: "OG", 1: 0})
    trepo["LGT"] = "trepo"
    df = pd.merge(df, trepo, how="inner").sort_values(by=["LGT", "OG", "Total"], ascending=[False, False, False])

    return df
