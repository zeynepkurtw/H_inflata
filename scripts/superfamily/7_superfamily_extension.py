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
    og_ann_aa_file = "data/superfamily/og_ann_aa.csv"

    out_file = "plots/family/heatmap_superfamily.png"


def plot_heatmap(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                           # "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[["C. membranifera", "K. bialata", "H. inflata", "S. salmonicida", "G. intestinalis",
             "G. muris"]]

    sns.set_theme(style="whitegrid")

    f, ax = plt.subplots(figsize=(10, 10))

    ax = sns.heatmap(df.T,
                     norm=LogNorm(),
                     cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     #annot=True,
                     cbar_kws={"shrink": .5}
                     ).set_title("superfamily_heatmap normalized & filtered")

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


# Common IPRs
def filter(df, filter):
    return df[df.apply(lambda c: (c[1:] > filter).any(), axis=1)].sort_values(by="HIN", ascending=False)


def normalization(df):
    df["HIN"] /= 142.6
    df["spiro"] /= 14.7
    df["wb"] /= 12.6
    df["muris"] /= 9.8
    df["carpe"] /= 24.2
    df["kbiala"] /= 51
    #trepo doesnt have a genome size
    return df


species = ['HIN', 'spiro', 'wb', 'muris', 'carpe', 'kbiala'] #trepo is excluded

df = pd.read_csv(og_ann_aa_file, header="infer", sep="\t")
df = df.dropna(subset="ipr").drop_duplicates(subset=["id", "ipr"])
count = dict_count_ipr(df, "HIN")

for sp in species:
    count = pd.merge(count, dict_count_ipr(df, sp), how="outer")

counts = count.set_index("ann_inter")
counts = filter(normalization(counts), 5 )
plot_heatmap( counts )
