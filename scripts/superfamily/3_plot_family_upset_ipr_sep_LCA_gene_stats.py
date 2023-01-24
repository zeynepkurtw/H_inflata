import pandas as pd
from upsetplot import UpSet
from upsetplot import plot
from matplotlib import pyplot as plt
import seaborn as sns
import glob

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    family_files = "data/superfamily/family_*"
    out_file = "plots/family/upset_family_{}_gene.png"

family = {}
list_files = glob.glob(family_files)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    i = i.split("_")[1]
    family[i] = element

def make_upset_data(df):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    #df = df.drop(columns=["db", "db_acc", "ann_db" ]).drop_duplicates()
    #df = df.drop_duplicates(subset = ["OG", "id", "ipr"])
    df = df.drop_duplicates(subset = ["id"])


    df_upset = df.set_index(df["C. membranifera"] >= 1). \
        set_index(df["K. bialata"] >= 1, append=True). \
        set_index(df["H. inflata"] >= 1, append=True). \
        set_index(df["Trepomonas pc1"] >= 1, append=True). \
        set_index(df["S. salmonicida"] >= 1, append=True). \
        set_index(df["G. intestinalis"] >= 1, append=True). \
        set_index(df["G. muris"] >= 1, append=True). \
        set_index(df["family"] == key, append=True)
    return df_upset

def upset_plot(df: pd.DataFrame, file_out: str, key: str) -> None:
    #sns.set_style("whitegrid", {'axes.grid': False})
    plt.figure(figsize=(10, 3))

    upset = UpSet(df,
                   intersection_plot_elements=20,
                   min_degree=2,
                   #min_subset_size=10,
                   show_counts=True,
                   sort_categories_by=None,
                   )

    upset.style_subsets(present=["family"],
                         facecolor="white",
                         edgecolor="black",
                         hatch="xxx",
                         linewidth=2,
                         label="family_{}_genes".format(key))

    upset.style_subsets(absent=["K. bialata", "C. membranifera"],
                         min_degree=5,
                         facecolor="blue",
                         label="Genes shared by diplomonad OGs")

    upset.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                         facecolor="red",
                         label="Genes shared by H. inflata and Trepomonas pc1 OGs")

    upset.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                         min_degree=4,
                         facecolor="green",
                         label="Genes shared by Free-living species OGs")

    upset.plot()
    #plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(file_out.format(key), format="png", bbox_inches='tight', dpi=600)
    plt.show()


def family_stats(df):
    len_fam_og = df.drop_duplicates("OG")
    len_fam_gene = df.drop_duplicates("id")
    len_fam_ipr = df
    len_fam_ipr_u = df.drop_duplicates("ipr")
    len_fam_pfam_u = df.groupby(["db"]).get_group( "Pfam").drop_duplicates("db_acc")

    print(key)
    print( "OG = " , len(len_fam_og))
    print( "genes = " , len(len_fam_gene))
    print( "IPR = " , len(len_fam_ipr))
    print( "IPR_unique = " , len(len_fam_ipr_u))
    print( "Pfam_unique = " , len(len_fam_pfam_u))
    print("")


sp_dict = {}
sp_list =  ['carpe', 'kbiala', 'HIN', 'trepo', 'spiro', 'wb', 'muris']

for key, value in family.items():
    family[key] = pd.read_csv(value, header="infer", sep="\t")
    df_upset = make_upset_data(family[key]).sort_values(by = ["OG", "id"])
    df_plot = upset_plot(df_upset, out_file, key)

    #family_stats(df_upset)
    for sp in sp_list:
        df = family[key].groupby("sp").get_group(sp)
        print(sp)
        family_stats(df)











