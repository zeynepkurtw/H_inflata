import pandas as pd
import glob
import seaborn as sns
import matplotlib
#matplotlib.use('Agg')  # Or try 'TkAgg', 'Qt5Agg', etc.
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_aa_file = "data/superfamily/4_og_ann_aa/og_ann_aa.csv"

    #out_file = "plots/family/4_heatmap_superfamily/heatmap_superfamily_IPR_gene_numbers.tiff"
    out_file = "plots/family/4_heatmap_superfamily/heatmap_superfamily_IPR_gene_numbers.png"


def plot_heatmap(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[["C. membranifera", "K. bialata", "H. inflata", "Trepomonas pc1", "S. salmonicida", "G. intestinalis",
             "G. muris"]]

    custom_labels = ["LRR", "HB", "CP", "PK", "CRP", "WD40", "ankyrin"]

    sns.set_theme(style="whitegrid")

    f, ax = plt.subplots(figsize=(10, 10))

    heatmap = sns.heatmap(df.T,
                     norm=LogNorm(),
                     cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     annot=True,
                     cbar_kws={"shrink": .5}
                     )
    #heatmap.set_title("Top7 Superfamily with Gene Numbers", y=1.1)
    heatmap.set_xlabel('')  # Removes the x-axis title

    # Positioning custom labels above each column
    for idx, label in enumerate(custom_labels):
        ax.text(x=idx + 0.5, y=-0.1, s=label, ha='center', va='center')

    #plt.savefig(out_file, format="tiff", bbox_inches='tight', dpi=1200)
    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
    plt.show()


def dict_count_ipr(df, sp) -> pd.DataFrame:
    """
    Count iprs for superfamily in each sp
    @return: superfamily count table
    """
    count = df.groupby(["sp"]).get_group(sp)
    count = count[count["ann_inter"].str.contains("superfamily", case=False)]
    count = count.groupby(["ipr", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)
    return count.rename(columns={0: sp})

species = ['HIN', 'trepo', 'spiro', 'wb', 'muris', 'carpe', 'kbiala'] #trepo is excluded

df = pd.read_csv(og_ann_aa_file, header="infer", sep="\t")
df = df.dropna(subset="ipr").drop_duplicates(subset=["id", "ipr"])
count = dict_count_ipr(df, "HIN")

for sp in species:
    count = pd.merge(count, dict_count_ipr(df, sp), how="outer")

counts = count.set_index("ipr")[:7]
plot_heatmap( counts )
