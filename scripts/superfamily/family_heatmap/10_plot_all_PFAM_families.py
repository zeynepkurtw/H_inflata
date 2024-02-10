import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


og_ann_ipr_file = "data/superfamily/5_og_ann_ipr/og_ann_ipr.csv"
out_file = "plots/family/6_heatmap_family_Pfam/heatmap_all_families_Pfam_entry.png"

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

    sns.set_style("whitegrid", {'axes.grid': False})

    plt.figure(figsize=(10, 30))

    # Create the heatmap
    ax = sns.heatmap(df,
                     norm=LogNorm(),
                     #cmap=sns.color_palette("light:b", as_cmap=True),
                     cmap= sns.color_palette("gray_r", as_cmap=True),
                     #cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=0.3,
                     linecolor="whitesmoke",
                     annot=True,
                     cbar_kws={"shrink": .3},
                     annot_kws={"size": 13})  # Adjust font size of annotations here

    ax.tick_params(axis='both', which='major', labelsize=15)  # Adjust the size as needed
    #remove x and y labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    y_labels = ax.get_yticklabels()
    ax.set_yticklabels(y_labels, rotation=0)
    ax.set_title("Pfam domain distribution for top 7 protein families")
    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
    return plt.show()
def plot_clustermap(df):

    sns.set_style("whitegrid", {'axes.grid': False})

    plt.figure(figsize=(10, 20))
    ax = sns.clustermap(df,
                    metric="correlation",
                    method="single",
                    standard_scale=1,
                     norm=LogNorm(),
                     cmap= sns.color_palette("Blues", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.2,
                     annot=True,
                     cbar_kws={"shrink": .5},
                     annot_kws={"size": 8})  # Adjust font size of annotations here

    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
    return plt.show()

def dict_count_ipr(df, sp) -> pd.DataFrame:
    """
    Count iprs for Pfam in each sp
    @return: PFAM count table
    """
    count = df.groupby(["sp"]).get_group(sp)
    count = count[count["db"].str.contains("Pfam", case=False)]
    count = count.groupby(["db_acc", "family"]).size().reset_index().sort_values(by=[0], ascending=False)
    return count.rename(columns={0: sp})

species = [ 'carpe', 'kbiala', 'HIN', 'trepo', 'spiro', 'wb', 'muris'] #trepo is excluded

df = pd.read_csv(og_ann_ipr_file, header="infer", sep="\t")
df = df.dropna(subset="db_acc").drop_duplicates(subset=["id", "db_acc"])
count = dict_count_ipr(df, "carpe")

for sp in species:
    count = pd.merge(count, dict_count_ipr(df, sp), how="outer")

counts = count.set_index("db_acc")
plot_heatmap( counts )















"""
# count ipr for each family
df = pd.read_csv(og_ann_ipr_file, header="infer", sep="\t")
df_pfam = df[df["db"] == "Pfam"]
df_pfam = df_pfam.groupby(["sp", "family", "ipr"]).size().reset_index().sort_values(by=["sp", "family", 0], ascending=False)
df_counts = df_pfam.rename(columns={0: "count"})
#df_counts['sp_family'] = df_counts['sp'] + '_' + df_counts['family']

df_pivot = df_counts.pivot(index="ipr", columns="sp", values="count").fillna(0)

plot_heatmap(df_pivot)
#plot_clustermap(df_pivot)
"""