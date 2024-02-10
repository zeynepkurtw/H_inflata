import pandas as pd
from upsetplot import UpSet
from matplotlib import pyplot as plt
import seaborn as sns
from utils import get_og_stack, get_og_count_sing

"""
OG gene list stack

plot only LGT OG and get genes from OG : the hashed bars in the previous version
"""

og_ann_lgt_file = "data/orthogroups/og_ann_lgt.csv"
og_ann_file = "data/orthogroups/og_ann.csv"

og_lgt_genes_csv = "data/orthogroups/og_lgt_genes.csv"
og_lgt_genes_excel = "data/orthogroups/og_lgt_genes.xlsx"

out = "plots/upset_trepo_LGT/upset_trepo_only_LGT_genes.png"


def make_upset_data(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1-LGT",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})
    # df = df.drop(columns= "Trepomonas pc1")
    df = df.drop_duplicates(subset=["id"])
    df = df.dropna(subset="LGT")

    # Number of cases where below conditions are true is plotted..
    # All trepo-LGT genes are marked as LGT
    # OGs that have these trepo-LGT genes are filtered
    # Genes that are found in these OGs together with the Trepo-LGT genes are shown in the upset as white
    # Genes that make up a specific sp combintions are colored as red, green, blue.

    df_upset = df.set_index(df["C. membranifera"] >= 1). \
        set_index(df["K. bialata"] >= 1, append=True). \
        set_index(df["H. inflata"] >= 1, append=True). \
        set_index(df["Trepomonas pc1-LGT"] >= 1, append=True). \
        set_index(df["S. salmonicida"] >= 1, append=True). \
        set_index(df["G. intestinalis"] >= 1, append=True). \
        set_index(df["G. muris"] >= 1, append=True)
    # set_index(df["LGT"] == "trepo", append=True) # read all the lines in a given df_upset with trepo in this case all the lines

    return df_upset


def upset_plot(file_out: str, df: pd.DataFrame) -> None:
    sns.set_style("whitegrid", {'axes.grid': False})
    fig = plt.figure(figsize=(10, 3))
    upset = UpSet(df,
                  intersection_plot_elements=10,
                  # element_size=60,
                  # fig=fig, element_size=None,
                  # min_degree=2,
                  # min_subset_size=5,
                  show_counts=True,
                  # sort_categories_by=None,
                  sort_categories_by="cardinality",
                  facecolor="darkgray",
                  sort_by="cardinality",
                  )

    upset.style_subsets(present=["Trepomonas pc1-LGT"],
                        facecolor="lightgray",

                        # hatch="xxx",
                        linewidth=2,
                        label="Metamonada")

    upset.style_subsets(absent=["K. bialata", "C. membranifera"],
                        min_degree=5,
                        facecolor="salmon",
                        label="Diplomonads")

    upset.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                        facecolor="darkseagreen",
                        label="2ndary Free-livings")

    upset.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                        min_degree=4,
                        facecolor="steelblue",
                        label="Free-livings")

    upset.plot()
    plt.suptitle("Intersections of LGT-Containing Orthogroups in Metamonada")
    plt.savefig(file_out, format="png", dpi=600)  # bbox_inches='tight', dpi=1200)
    plt.show()


def get_og_ann_contains_lgt(df, og_ann):
    df = df["OG"].drop_duplicates()
    return pd.merge(og_ann, df, on="OG", how="inner")


# Read og_ann_lgt
og_ann_lgt = pd.read_csv(og_ann_lgt_file, header="infer", sep="\t").dropna(subset="LGT")
og_ann = pd.read_csv(og_ann_file, header="infer", sep="\t")

# Make upset df and rename columns
df_upset = make_upset_data(og_ann_lgt)

og_ann_lgt_genes = get_og_ann_contains_lgt(df_upset, og_ann)
og_ann_lgt_genes.to_csv(og_lgt_genes_csv, sep="\t", index=False)
og_ann_lgt_genes.to_excel(og_lgt_genes_excel, index=False)

# Plot upset
df_plot = upset_plot(out, df_upset)

print("# of OG in df_upset", len(df_upset["OG"].drop_duplicates()))
print("# of genes in df_upset", len(df_upset["id"].drop_duplicates()))
print(" ")
print("# of OG in og_ann_lgt_genes", len(og_ann_lgt_genes["OG"].drop_duplicates()))
print("# of genes in og_ann_lgt_genes", len(og_ann_lgt_genes["id"].drop_duplicates()))

og_ann_lgt_genes = og_ann_lgt_genes.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1-LGT",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

trepo= og_ann_lgt_genes[(og_ann_lgt_genes['Trepomonas pc1-LGT'] == True)].drop_duplicates(subset= "id")
hin= og_ann_lgt_genes[(og_ann_lgt_genes['H. inflata'] == True)].drop_duplicates(subset= "id")
spiro= og_ann_lgt_genes[(og_ann_lgt_genes['S. salmonicida'] == True)].drop_duplicates(subset= "id")
wb= og_ann_lgt_genes[(og_ann_lgt_genes['G. intestinalis'] == True)].drop_duplicates(subset= "id")
muris= og_ann_lgt_genes[(og_ann_lgt_genes['G. muris'] == True)].drop_duplicates(subset= "id")

