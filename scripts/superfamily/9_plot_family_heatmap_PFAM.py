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

    out_file = "plots/family/heatmap_family_{}_Pfam.png"


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

    f, ax = plt.subplots(figsize=(20, 10))

    ax = sns.heatmap(df.T,
                     norm=LogNorm(),
                     cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     annot=True,
                     cbar_kws={"shrink": .5}
                     ).set_title("family_{}_heatmap".format(family))

    plt.savefig(out_file.format(family), format="png", bbox_inches='tight', dpi=1200)
    return plt.show()


def dict_count_ipr(df, family, sp) -> pd.DataFrame:
    """
    Count iprs for superfamily in each sp
    @param df: og_ann_ipr lecuine, cysteine, ankyrin
    @param fam: family name : str
    @param sp: species : str
    @return: family dict with ipr, ann_ipr and counts
    """
    count = df.groupby(["family", "sp"]).get_group((family, sp))
    count = count[count["db"] == "Pfam"]
    count = count.groupby(["db_acc", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)
    return count.rename(columns={0: sp})


families = {
    "fam_l": "leucine rich repeat",
    "fam_c": "cysteine rich",
    "fam_a": "ankyrin repeat"
}
species = ['HIN','trepo', 'spiro', 'wb', 'muris', 'carpe', 'kbiala']
counts = {}

df = pd.read_csv(og_ann_ipr_LCA_file, header="infer", sep="\t")
for fam, aa in families.items():
    count = dict_count_ipr(df, aa, "HIN")

    for sp in species:
        count = pd.merge(count, dict_count_ipr(df, aa, sp), how="outer")

    counts[fam] = count.set_index(["db_acc","ann_inter"])
    plot_heatmap(counts[fam], fam)
