import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Define paths for CD-hit fasta files and InterPro annotation files
path_cdhit = "/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/cdhit_fasta/*"
path_interpro = "/Users/zeyku390/PycharmProjects/H_inflata/data/interpro_ann/*.csv"

out1 = "/Users/zeyku390/PycharmProjects/H_inflata/plots/Figure1/figure1b_functional_h.png"
out2 = "/Users/zeyku390/PycharmProjects/H_inflata/plots/Figure1/figure1b_unique_h.png"

# Dictionary with total genes per species
dic_total_genes = {
    "HIN": 79341,
    "spiro": 8661,
    "wb": 4963,
    "muris": 4653,
    "carpe": 8300,
    "kbiala": 17389,
    "trepo": 7980
}

# Function to extract IDs from fasta files
def extract_fasta_ids(fasta_path):
    ids = [record.id for record in SeqIO.parse(fasta_path, "fasta")]
    return pd.DataFrame(ids, columns=['id'])

# Function to calculate value counts from a dataframe
def calculate_value_counts(df):
    melted_df = df.melt(var_name='columns', value_name='values').dropna().drop_duplicates()
    count_df = melted_df['columns'].value_counts().reset_index()
    count_df.columns = ['columns', 'counts']
    return count_df

# Function to read CD-hit files and InterPro annotation files
def read_files_cdhit(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split("/")[-1]
        i = i.split(".cdhit")[0]
        dic[i] = element
    return dic

def read_files_interpro(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split("/")[-1]
        i = i.split(".csv")[0]
        dic[i] = element
    return dic

# Main function to get functional CD-hit clusters
def get_functional_cdhit_clusters(dic_cd, dic_int):
    id_dict = {}
    annotation_dict = {}

    # Iterate over CD-hit files and corresponding InterPro files
    for identifier, cdhit_path in dic_cd.items():
        id_dict[identifier] = extract_fasta_ids(cdhit_path)
        for sp, interpro_path in dic_int.items():
            if sp in identifier:
                df_int = pd.read_csv(dic_int[sp], sep="\t", header="infer").drop_duplicates("id")["id"]
                annotation_dict[identifier] = pd.merge(id_dict[identifier], df_int, on="id", how="inner")

    functional_dict = annotation_dict.copy()

    # Concatenate dataframes
    df_id = pd.concat(id_dict, axis=1)
    df_functional = pd.concat(functional_dict, axis=1)

    # Calculate counts for unique and functional groups
    unique_groups = calculate_value_counts(df_id).rename(columns={"counts": "unique"})
    functional_groups = calculate_value_counts(df_functional).rename(columns={"counts": "functional"})

    # Prepare data for plotting
    functional_groups = functional_groups["columns"].str.split("_", expand=True).rename(columns={0: "sp", 1: "cdhit"}).join(functional_groups["functional"])
    df_plot = functional_groups.join(unique_groups["unique"])

    species_mapping = {
        "carpe": "C. membranifera",
        "kbiala": "K. bialata",
        "HIN": "H. inflata",
        "trepo": "Trepomonas pc1",
        "spiro": "S. salmonicida",
        "wb": "G. intestinalis",
        "muris": "G. muris"
    }

    df_plot['subtract'] = df_plot.apply(lambda row: dic_total_genes[row['sp']] - row['unique'], axis=1)

    return df_plot

# Load CD-hit and InterPro files
dic_cd = read_files_cdhit(path_cdhit)
dic_int = read_files_interpro(path_interpro)

# Generate plot data
df_plot = get_functional_cdhit_clusters(dic_cd, dic_int)

# Pivot the data for heatmap
df_heatmap = df_plot.pivot(index="sp", columns="cdhit", values="functional")

# Create and display the heatmap
sns.set_style("whitegrid", {'axes.grid': False})
plt.figure(figsize=(10, 8))

ax = sns.heatmap(df_heatmap,
                 cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                 annot=True, fmt="d")
ax.set_xlabel("Sequence Identity")
ax.set_ylabel("Species")
ax.set_title("Heatmap of Functional Clusters")

plt.savefig(out1, format="png", dpi=1200)
plt.show()
