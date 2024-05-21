import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns

"""
Visualization of Cd-hit clustering
Input: cdhit.fasta
Output: df_func, df_hypo, df_id, df_plot, and a line plot
"""

# Load files from specified directory
# Uncomment the relevant path as needed
#path_cdhit = '/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/*.cdhit.clstr'
#path_cdhit = '/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/cdhit_slow/*.cdhit'
path_cdhit = "/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/cdhit_fasta/*"
path_interpro = "/Users/zeyku390/PycharmProjects/H_inflata/data/interpro_ann/concat/ann_ipr_cat.csv"

def load_file_paths(pattern):
    return {file.split("/")[-1].split(".cdhit")[0]: file for file in glob.glob(pattern)}

def extract_ids_from_fasta(fasta_file):
    ids = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    return pd.DataFrame(ids, columns=['id'])

# Merge IDs with InterPro annotations
file_paths = load_file_paths(path_cdhit)
df_annotations = pd.read_csv(path_interpro, sep="\t", header="infer").drop_duplicates("id")

#merged_data = {key: pd.merge(extract_ids_from_fasta(path_cdhit), df_annotations, how="left") for key, path in file_paths.items()}

def merge_ids_with_annotations(file_paths, annotations_df):
    """
    Merges IDs extracted from fasta files with annotations.

    :param file_paths: A dictionary of file paths to fasta files, where each key is a unique identifier.
    :param annotations_df: A DataFrame containing annotations to be merged with the extracted IDs.
    :return: A dictionary of DataFrames, each representing merged data for a file.
    """
    merged_data = {}
    for key, fasta_path in file_paths.items():
        # Extract IDs from the current fasta file
        ids_df = extract_ids_from_fasta(fasta_path)

        # Merge the extracted IDs with the annotations DataFrame
        merged_df = pd.merge(ids_df, annotations_df, how="left", on="id")

        # Store the merged DataFrame in the dictionary using the key
        merged_data[key] = merged_df

    return merged_data


# Now, you can call this function using the file_paths dictionary and the annotations DataFrame
merged_data = merge_ids_with_annotations(file_paths, df_annotations)


def perform_anti_join(x, y):
    """Exclude rows in x that are present in y"""
    return pd.merge(x, y, how='left', indicator=True).query('_merge == "left_only"').drop(columns='_merge')

# Distinguish between functional and hypothetical proteins
functional_data, hypothetical_data = {}, {}
for key, df in merged_data.items():
    functional_data[key] = df
    #hypothetical_data[key] = perform_anti_join(extract_ids_from_fasta(file_paths[key]), df)

# Consolidate dictionaries into DataFrames
#df_id = pd.concat(functional_data, axis=1).droplevel(1, axis=1)
#df_hypo = pd.concat(hypothetical_data, axis=1).droplevel(1, axis=1)
df_func = pd.concat(functional_data, axis=1).droplevel(1, axis=1)

def count_values(df):
    counts = df.melt(var_name='Type', value_name='Counts').dropna().drop_duplicates()['Type'].value_counts().reset_index()
    return counts

# Prepare data for plotting
#id_counts = count_values(df_id).rename(columns={"index": "Identity", "Type": "ID"})
func_counts = count_values(df_func).rename(columns={"index": "Identity", "Type": "Functional"})
# Data for plotting
#plot_data = pd.concat([id_counts, func_counts], axis=0)

# Create and customize the plot
sns.set_style("whitegrid", {'axes.grid': False})
palette = sns.dark_palette("#69d", reverse=False, n_colors=2)
#figure = sns.lineplot(data=plot_data, x="Identity", y="Counts", hue="Type", palette=palette)
figure = sns.lineplot(data=func_counts, x="Identity", y="Counts", hue="Type", palette=palette)

figure.set_xticklabels(['70%', '75%', '80%', '85%', '90%', '95%', '100%'])
figure.set_xlabel("Sequence Identity")
figure.set_ylabel("Number of Clusters")
figure.legend(title='Cluster Type', loc='upper left', labels=['Unique', 'Functional'])

# Display or save the plot
line_plot = figure.get_figure()
#line_plot.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1b.png', dpi=1200)
line_plot.show()