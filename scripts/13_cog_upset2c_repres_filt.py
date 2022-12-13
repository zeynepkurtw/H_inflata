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
    og_count = "data/orthogroups/og_count.csv"
    og_ann = "data/orthogroups/og_ann.csv"

    out_file = "plots/upset_cog2c_repres_filt.png"


def make_upset_data(og_count, og_ann):
    #df = df.drop(columns= ["OG", "Total"])
    #df.at["type"] = ["a", "b"]
    #df.loc[:, 'type'] = "type"
    og_ann = og_ann.dropna(subset="COG")
    og_ann["COG"] = og_ann["COG"].str.split(",")
    og_ann = og_ann.explode("COG").reset_index(drop=True)
    og_ann = og_ann[(og_ann["COG"] != "S")]
    og_ann = og_ann[(og_ann["COG"] != "R")]

    # Count number of genes that belongs to a COG categ. in each OG
    # Take the most abundant COG as a representative of teh group
    # UPset plot return counts for number of genes on the stacked bars.
    df_cog =  og_ann.groupby(["OG", "COG"]).agg(count=('COG', 'size')).\
        reset_index().sort_values(by=["OG", "count"], ascending=[True, False])
    df= pd.merge(df_cog, og_count)
    df= df.drop_duplicates(subset="OG", keep="first")

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
        set_index(df["G. muris"] >= 1, append=True)
    # set_index(df["LGT"] == "trepo", append=True)
    return df_upset, df_cog


def upset_plot(file_out: str, df: pd.DataFrame) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    fig = plt.figure(figsize=(10, 3))
    upset = UpSet(df,
                  intersection_plot_elements=0,
                  min_degree=2,
                  max_degree=5,
                  min_subset_size=20,
                  show_counts=True,
                  sort_categories_by=None)

    #pal = sns.color_palette("Set3")
    pal = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
           '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99',
           ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69',
            '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']]
    upset.add_stacked_bars(by="COG",
                           sum_over="count",
                           colors=pal,
                           title="Total number of OG for each COG cat.",
                           elements=10)

    upset.style_subsets(absent=["K. bialata", "C. membranifera"],
                        min_degree=5,
                        facecolor="blue",
                        label="OGs shared by all diplomonads")


    upset.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                        min_degree=2,
                        facecolor="green",
                        label="OGs shared by Free-living species")

    upset.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                        facecolor="red",
                        label="OGs shared by H. inflata and Trepomonas pc1")

    upset.plot()

    plt.suptitle("COG categories og>20")

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.savefig(out_file, format="png", dpi=1200, bbox_inches='tight')
    plt.show()



og_ann= pd.read_csv(og_ann, sep="\t", header='infer',
                 dtype={"carpe": int, "kbiala": int, "HIN": int, "trepo": int, "wb": int, "spiro": int, "muris": int,
                        "Total": int, "OG": str, "id": str, "ann_f": str, "ann_inter": str, "ann_egg": str, "ipr": str,
                        "COG": str, "KEGG_KOs": str, "LGT": str, "Type": str})

og_count = pd.read_csv(og_count, sep="\t", header='infer')


# Make upset df and rename columns
df_upset, df_cog = make_upset_data(og_count, og_ann)

# Plot upset
df_plot = upset_plot(out_file, df_upset)
