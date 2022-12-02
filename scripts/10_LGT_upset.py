import pandas as pd
from upsetplot import UpSet
from matplotlib import pyplot as plt
import seaborn as sns
from utils import get_og_stack, get_og_count_sing

"""
OG gene list stack
"""

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    trepo_lgt_file = 'resource/5_LGT/trepo_list.xlsx'
    out = "plots/upset_trepo.png"

"""
Find Trepo LGT list in OG
"""
def trepo_og():
    trepo = pd.read_excel(trepo_lgt_file, header=None, index_col=0, skiprows=3).reset_index()
    df = pd.merge(trepo, df_og, on=0)[[0,1,"OG"]]
    df["LGT"] = "trepo"
    return df


def trepo_upset(df_trepo_og, df_count_s):
    df = pd.merge(df_trepo_og, df_count_s, on="OG",  how="outer")#.drop(columns="OG")
    return df


def upset_data(df_count_s):
    df_count_s = df_count_s.rename(columns={"carpe": "C. membranifera",
                                            "kbiala": "K. bialata",
                                            "HIN": "H. inflata",
                                            "trepo": "Trepomonas pc1",
                                            "spiro": "S. salmonicida",
                                            "wb": "G. intestinalis",
                                            "muris": "G. muris"})

    df_upset = df_count_s.set_index(df_count_s["C. membranifera"] >= 1). \
        set_index(df_count_s["K. bialata"] >= 1, append=True). \
        set_index(df_count_s["H. inflata"] >= 1, append=True). \
        set_index(df_count_s["Trepomonas pc1"] >= 1, append=True). \
        set_index(df_count_s["S. salmonicida"] >= 1, append=True). \
        set_index(df_count_s["G. intestinalis"] >= 1, append=True). \
        set_index(df_count_s["G. muris"] >= 1, append=True). \
        set_index(df_count_s["LGT"] == "trepo", append=True)
    return df_upset


def upset_plot(file_out: str, df: pd.DataFrame) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    upset3 = UpSet(df,
                   intersection_plot_elements=10,
                   min_degree=2,
                   # min_subset_size=10,
                   show_counts=True,
                   # sort_categories_by="cardinality",
                   # sort_by="cardinality",
                   )

    upset3.style_subsets(present=["LGT"],
                         facecolor="pink",
                         label="trepo putative LGT genes")

    upset3.style_subsets(present=["LGT"],
                         min_subset_size=10,
                         label="> 10", hatch="xx",
                         edgecolor="black")

    upset3.plot()
    # plt.suptitle("OGs shared by all diplomonads")
    plt.savefig(file_out, format="png",
                dpi=600)  # bbox_inches='tight', dpi=1200)
    plt.show()


df_og = get_og_stack()
df_trepo_og = trepo_og()
df_count_s = get_og_count_sing()
df_upset = trepo_upset(df_trepo_og, df_count_s)
df_plot = upset_data(df_upset)
upset_plot(out, df_plot)