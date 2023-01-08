import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_ipr_LCA_file = "data/superfamily/og_ann_ipr_LCA.csv"

    out_file = "plots/family/heatmap_family_{}.png"


def plot_heatmap(df, family):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[["C. membranifera", "K. bialata", "H. inflata", "Trepomonas pc1", "S. salmonicida", "G. intestinalis",
             "G. muris"]]

    sns.set_theme(style="whitegrid")

    f, ax = plt.subplots(figsize=(10, 10))

    plot = sns.heatmap(df.T,
                       norm=LogNorm(),
                       cmap=sns.color_palette("ch:start=.2,rot=-.3",
                                              as_cmap=True),
                       square=True,
                       cbar=False,
                       fmt='g',
                       linewidths=.4,
                       cbar_kws={"shrink": .5}).set_title(
        "family_{}_heatmap".format(family))

    plt.savefig(out_file.format(key), format="png", dpi=600)  # bbox_inches='tight', dpi=1200)
    return plt.show()


def dict_count_ipr(df, fam_dict, fam, sp):
    """
    Count iprs for superfamily in each sp
    @param df: og_ann_ipr lecuine, cysteine, ankyrin
    @param fam_dict: empty family dictionary
    @param fam: family name : str
    @param sp: species : str
    @return: family dict with ipr, ann_ipr and counts
    """
    fam_dict={}
    fam_dict[sp] = df.groupby(["family", "sp"]).get_group((fam, sp))
    fam_dict[sp] = fam_dict[sp].groupby(["ipr", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)
    fam_dict[sp] = fam_dict[sp].rename(columns={0: sp})

    return fam_dict


df = pd.read_csv(og_ann_ipr_LCA_file, header="infer", sep="\t")


fam = {
    "fam_l": "leucine rich repeat",
    "fam_c": "cysteine rich",
    "fam_a": "ankyrin repeat"
}
sp_list = ['HIN', 'trepo', 'spiro', 'wb', 'muris', 'carpe', 'kbiala']

for sp in sp_list:
    for family, aa in fam.items():
        fam_dict = dict_count_ipr(df, family, aa, sp)

        for sp, count in fam_dict.items():
            df = fam_dict["HIN"]

            if sp != "HIN":
                df = pd.merge(df, count, how="outer").set_index("ann_inter")

        plot_heatmap(df, aa)

"""
family_list = ["cysteine rich", "leucine rich repeat", "ankyrin repeat"]
sp_list = ['HIN', 'trepo', 'spiro', 'wb', 'muris', 'carpe', 'kbiala' ]
family = {}
all_d = {}
pfam = {}
fam_c = {}
fam_l = {}
fam_a = {}

df = pd.read_csv(og_ann_ipr_LCA_file, header="infer", sep="\t")
for i in family_list:
    for j in sp_list:
        family[i,j] = df.groupby(["family", "sp"]).get_group((i, j))
        all_d[i,j] = family[i,j].groupby(["ipr", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)

        pfam[i,j] = family[i,j][family[i,j]["db"] == "Pfam"]
        pfam[i,j] = pfam[i,j].groupby(["ipr", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)
"""

# plot
