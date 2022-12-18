import pandas as pd
from upsetplot import UpSet, plot
from matplotlib import pyplot as plt
import seaborn as sns
from utils import *

"""
OG gene list stack
"""

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    trepo_lgt_file = 'resource/5_LGT/trepo_list.xlsx'
    out = "plots/upset_trepo_stacked.png"


def merge_trepo_og(trepo_lgt_file, df_og, df_count_s):
    trepo = pd.read_excel(trepo_lgt_file, header=None, index_col=0, skiprows=3).reset_index()
    df = pd.merge(trepo, df_og, on=0)
    df["LGT"] = "LGT-Trepomonas pc1"
    df = df[["OG", 0, 1, "LGT"]]

    df = pd.merge(df, df_count_s, on="OG", how="outer")
    return df


def make_upset_data(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df_upset = df.set_index(df["C. membranifera"] >= 1). \
        set_index(df["K. bialata"] >= 1, append=True). \
        set_index(df["H. inflata"] >= 1, append=True). \
        set_index(df["Trepomonas pc1"] >= 1, append=True). \
        set_index(df["S. salmonicida"] >= 1, append=True). \
        set_index(df["G. intestinalis"] >= 1, append=True). \
        set_index(df["G. muris"] >= 1, append=True)  # . \
    # set_index(df["LGT"] == "trepo", append=True)

    df_upset["LGT"] = df_upset["LGT"].fillna("")

    return df_upset


def upset_plot(file_out: str, df: pd.DataFrame) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    upset = UpSet(df,
                  # intersection_plot_elements=10,
                  intersection_plot_elements=0,
                  min_degree=2,
                  min_subset_size=10,
                  show_counts=True,
                  #element_size=60,
                  #intersection_plot_elements=10,
                  sort_categories_by=None,
                  # sort_categories_by="cardinality",
                  # sort_by="cardinality",
                  )

    pal = sns.dark_palette("#69d", n_colors=2, reverse=False)
    upset.add_stacked_bars(by="LGT",
                           colors=["black", "pink"],
                           title="Intersection size",
                           elements=10)

    upset.style_subsets(absent=["K. bialata", "C. membranifera"],
                        min_degree=5,
                        facecolor="blue",
                        label="OGs shared by all diplomonads")

    upset.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                        facecolor="red",
                        label="OGs shared by H. inflata and Trepomonas pc1")

    upset.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                        min_degree=4,
                        facecolor="green",
                        label="OGs shared by Free-living species")
    upset.plot()

    # plt.suptitle("OGs shared by all diplomonads")
    plt.savefig(file_out, format="png", dpi=600)  # bbox_inches='tight', dpi=1200)
    plt.show()


# Get og and genes
df_og = get_og_stack()

# Get og count table together with singletons
# df_count_s = get_og_count_sing()

# Get og count table together without singletons
df_count_s = get_og_count()

# Add trepo lgt genes
df_trepo = merge_trepo_og(trepo_lgt_file, df_og, df_count_s)
df_trepo.to_csv("data/LGT/trepo_lgt.csv", sep="\t", index=False, header=None)

# Make upset df and rename columns
df_upset = make_upset_data(df_trepo)

# Plot upset
df_plot = upset_plot(out, df_upset)

# upset_plot(out, df_plot)
