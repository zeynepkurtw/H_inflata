import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_aa_file = "data/superfamily/4_og_ann_aa/og_ann_aa.csv"

    out_file = "data/superfamily/stats.csv"


def all_db_stats(df):
    df = df.drop(columns=["db", "db_acc", "ann_db"]).drop_duplicates()
    df = df.dropna(subset="ipr")

    len_og = df.drop_duplicates("OG")
    len_gene = df.drop_duplicates("id")
    len_ipr = df
    len_ipr_u = df.drop_duplicates("ipr")

    print("IPR statistics")
    print("OG = ", len(len_og))
    print("genes = ", len(len_gene))
    print("IPR = ", len(len_ipr))
    print("IPR_unique = ", len(len_ipr_u))
    return df


def Pfam_stats(df):
    df = df[df["db"] == "Pfam"]
    df = df.dropna(subset="db_acc")

    len_og = df.drop_duplicates("OG")
    len_gene = df.drop_duplicates("id")
    len_pfam = df
    len_pfam_u = df.drop_duplicates("db_acc")

    print("Pfam statistics")
    print("OG = ", len(len_og))
    print("genes = ", len(len_gene))
    print("Pfam = ", len(len_pfam))
    print("Pfam_unique = ", len(len_pfam_u))
    return df


sp_list = ['carpe', 'kbiala', 'HIN', 'trepo', 'spiro', 'wb', 'muris', ]

og_ann_aa = pd.read_csv(og_ann_aa_file, header="infer", sep="\t")
for sp in sp_list:
    df = og_ann_aa.groupby("sp").get_group(sp)

    print(sp)
    df_all = all_db_stats(df)
    print("")
    df_pfam = Pfam_stats(df)
    print("")

