import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


og_ann_ipr_file = "data/superfamily/5_og_ann_ipr/plot_pfam_heatmap/og_ann_ipr.csv"
out_file = "plots/family/6_heatmap_family_Pfam/plot_pfam_heatmap/heatmap_all_families_Pfam_entry_filt.png"

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

    plt.figure(figsize=(5, 15))

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
                     annot_kws={"size": 8})  # Adjust font size of annotations here

    ax.tick_params(axis='both', which='major', labelsize=8)  # Adjust the size as needed
    #remove x and y labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    y_labels = ax.get_yticklabels()
    ax.set_yticklabels(y_labels, rotation=0)
    #ax.set_title("Pfam domain distribution for top 7 protein families filtered")
    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=800)
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

def filter(df, filter):
    return df[df.apply(lambda c: (c[1:] > filter).any(), axis=1)].sort_values(by="family", ascending=False)

# make a function that filters teh top 5 pfam fro each family
def filter_top5(df, filter):
    return df.groupby("family").head(filter).sort_values(by=["HIN", "family"], ascending=False)

species = [ 'carpe', 'kbiala', 'HIN', 'trepo', 'spiro', 'wb', 'muris']

df = pd.read_csv(og_ann_ipr_file, header="infer", sep="\t")
df = df.dropna(subset="db_acc").drop_duplicates(subset=["id", "db_acc"])
count = dict_count_ipr(df, "carpe")

for sp in species:
    count = pd.merge(count, dict_count_ipr(df, sp), how="outer")

counts = count.set_index("db_acc")
#counts = filter( counts, 10 )# domains with at least 10 proteins
counts = filter_top5( counts,5 )# domains with at least 10 proteins

plot_heatmap( counts )













