import pandas as pd
from upsetplot import UpSet
from matplotlib import pyplot as plt
import seaborn as sns
import glob

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_aa_file = "data/superfamily/og_ann_aa.csv"
    out_file = "plots/family/upset_family_{}.png"

aa_percs = ["perc_L", "perc_C"]


def make_upset_data(df, aa_perc):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df.drop_duplicates(subset=["id"]).sort_values(by= ["OG", "id"])

    df_upset = df.set_index(df["C. membranifera"] >= 1). \
        set_index(df["K. bialata"] >= 1, append=True). \
        set_index(df["H. inflata"] >= 1, append=True). \
        set_index(df["Trepomonas pc1"] >= 1, append=True). \
        set_index(df["S. salmonicida"] >= 1, append=True). \
        set_index(df["G. intestinalis"] >= 1, append=True). \
        set_index(df["G. muris"] >= 1, append=True). \
        set_index(df[aa_perc] > threshold, append=True)
    return df_upset


def upset_plot(df: pd.DataFrame, file_out: str, aa_perc: str) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    fig = plt.figure(figsize=(10, 3))
    upset = UpSet(df,
                  intersection_plot_elements=20,
                  min_degree=2,
                  min_subset_size=50,
                  show_counts=True,
                  sort_categories_by=None,
                  )

    upset.style_subsets(present=[aa_perc],
                        facecolor="white",
                        edgecolor="black",
                        hatch="xxx",
                        linewidth=2,
                        label="family_{}_{} >50".format(aa_perc, threshold))

    upset.style_subsets(absent=["K. bialata", "C. membranifera"],
                        min_degree=5,
                        facecolor="blue",
                        label="Genes shared by all diplomonads")

    upset.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                        facecolor="red",
                        label="Genes shared by H. inflata and Trepomonas pc1")

    upset.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                        min_degree=4,
                        facecolor="green",
                        label="Genes shared by Free-living species")

    upset.plot()
    plt.savefig(file_out.format(aa_perc), format="png", dpi=900)  # bbox_inches='tight', dpi=1200)
    plt.show()

def stats(df, aa_perc, threshold):
    df = df.drop(columns=["db", "db_acc", "ann_db"]).drop_duplicates()
    df =df[df[aa_perc] > threshold ]
    len_fam_og = df.drop_duplicates("OG")
    len_fam_gene = df.drop_duplicates("id")
    len_fam_ipr = df
    len_fam_ipr_u = df.drop_duplicates("ipr")

    print(aa_perc)
    print("OG = ", len(len_fam_og))
    print("genes = ", len(len_fam_gene))
    print("IPR = ", len(len_fam_ipr))
    print("IPR_unique = ", len(len_fam_ipr_u))
    print("")




for aa_perc in aa_percs:
    df = pd.read_csv(og_ann_aa_file, header="infer", sep="\t")
    if aa_perc == "perc_L":
        threshold = 10
        df_upset = make_upset_data(df, aa_perc)
        df_plot = upset_plot(df_upset, out_file, aa_perc)
        stats(df, aa_perc, threshold)

    elif aa_perc == "perc_C":
        threshold = 5
        df_upset = make_upset_data(df, aa_perc)
        df_plot = upset_plot(df_upset, out_file, aa_perc)
        stats(df, aa_perc, threshold)


