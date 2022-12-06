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
    og_ann_lgt = "data/orthogroups/og_ann_lgt.csv"

    out_file = "plots/upset_cog_lgt.png"


def make_upset_data(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[df["LGT"] == "trepo"]
    df = df.dropna(subset=['COG'])
    df["COG"] = df["COG"].apply(lambda x: x.split(",")[0])
    df = df.drop_duplicates(subset=["id"])

    df_upset = df.set_index(df["C. membranifera"] >= 1). \
        set_index(df["K. bialata"] >= 1, append=True). \
        set_index(df["H. inflata"] >= 1, append=True). \
        set_index(df["Trepomonas pc1"] >= 1, append=True). \
        set_index(df["S. salmonicida"] >= 1, append=True). \
        set_index(df["G. intestinalis"] >= 1, append=True). \
        set_index(df["G. muris"] >= 1, append=True)
    # set_index(df["LGT"] == "trepo", append=True)
    return df_upset


def upset_plot(file_out: str, df: pd.DataFrame) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    fig = plt.figure(figsize=(10, 3))
    upset = UpSet(df,
                  intersection_plot_elements=0,
                  min_degree=1,
                  #min_subset_size=10,
                  show_counts=True,
                  sort_categories_by=None)

    # pal = sns.dark_palette("#69d", n_colors=2, reverse=False)
    # pal = sns.dark_palette("#69d", reverse=True)
    pal = sns.color_palette("Set3")
    # pal=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']

    upset.add_stacked_bars(by="COG",
                           colors=pal,
                           title="Intersection size",
                           elements=10)

    upset.plot()

    # plt.suptitle("OGs shared by all diplomonads")

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.savefig(out_file, format="png", dpi=1200, bbox_inches='tight')
    plt.show()


df = pd.read_csv(og_ann_lgt, sep="\t", header='infer',
                 dtype={"carpe": int, "kbiala": int, "HIN": int, "trepo": int, "wb": int, "spiro": int, "muris": int,
                        "Total": int, "OG": str, "id": str, "ann_f": str, "ann_inter": str, "ann_egg": str, "ipr": str,
                        "COG": str, "KEGG_KOs": str, "LGT": str, "Type": str})

# Make upset df and rename columns
df_upset = make_upset_data(df)

# Plot upset
df_plot = upset_plot(out_file, df_upset)
