import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns
import matplotlib as plt

"""
Cdhit line plot
input= cdhit.fasta
output= df_func, df_hypo, df_id, df_plot, line_plot
"""

path_cdhit = "/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/cdhit_fasta/*"
#path_interpro= "/Users/zeyku390/PycharmProjects/H_inflata/data/interpro_ann/concat/ann_ipr_cat.csv"
path_interpro= "/Users/zeyku390/PycharmProjects/H_inflata/data/interpro_ann/*.csv"


def get_id(fasta_path):
    id= []
    with open(fasta_path, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
    return pd.DataFrame(id)

"get value counts from dataframe"
def get_counts(df):
    df_melt=df.melt( var_name='cols', value_name='vals').dropna().drop_duplicates()
    df_count=df_melt["cols"].value_counts().reset_index()
    return df_count

def get_functional_cdhit_clusters_(cdhit_files, interpro_file):
    dic_id={}
    dic_ann={}

    df_ann= pd.read_csv(path_interpro, sep="\t", header=None).drop_duplicates(0) #get unique genes
    for key, values in read_files_cdhit(cdhit_files).items():
        dic_id[key] = get_id(values)
        dic_ann[key]= pd.merge(dic_id[key], df_ann)

    dic_func = {}
    for key, values in dic_ann.items():
        dic_func[key] = values

    "concatanate dic elements in to a dataframe"
    df_id=pd.concat(dic_id, axis=1)
    df_id.columns= df_id.columns.droplevel(1) #remove multindex columns

    df_func=pd.concat(dic_func, axis=1)
    df_func.columns= df_func.columns.droplevel(1) #remove multindex columns

    id_groups = get_counts(df_id).rename(columns={"cols": "unique"})
    func_groups = get_counts(df_func).rename(columns={"cols": "functional"})
    #df_plot = pd.concat([id_groups, func_groups], axis=0)
    #df_plot_melt = df_plot.melt("index", var_name='cols', value_name='vals').dropna().sort_values(by="index",ascending=True)

    """ plot only the functional clusters"""
    df_plot_melt = func_groups.melt("index", var_name='cols', value_name='vals').dropna().sort_values(by="index",ascending=True)

    df_plot_melt["index"]= df_plot_melt["index"].apply( lambda x: x.split("_")[1])
    return df_plot_melt

def get_functional_cdhit_clusters(cdhit_files, interpro_files):
    """
    Merges CD-hit clusters with InterPro annotations for multiple species,
    extracting functional clusters.

    :param cdhit_files: Dictionary of CD-hit file paths keyed by species or identifier.
    :param interpro_files: Dictionary of InterPro annotation file paths keyed by the same species or identifier.
    :return: DataFrame of melted functional groups for plotting.
    """
    dic_id = {}
    dic_ann = {}

    # Loop through each species or identifier
    for key, cdhit_path in cdhit_files.items():
        # Load CD-hit IDs
        dic_id[key] = get_id(cdhit_path)

        # Load corresponding InterPro annotations for the species
        df_ann = pd.read_csv(interpro_files[key], sep="\t", header=None).drop_duplicates(0)

        # Merge CD-hit IDs with InterPro annotations
        dic_ann[key] = pd.merge(dic_id[key], df_ann, how="left")

    dic_func = {key: values for key, values in dic_ann.items()}

    # Concatenate dictionary elements into a DataFrame
    df_id = pd.concat(dic_id, axis=1)
    df_id.columns = df_id.columns.droplevel(1)  # Remove multi-index columns

    df_func = pd.concat(dic_func, axis=1)
    df_func.columns = df_func.columns.droplevel(1)  # Remove multi-index columns

    # Calculate counts for unique and functional groups
    id_groups = get_counts(df_id).rename(columns={"cols": "unique"})
    func_groups = get_counts(df_func).rename(columns={"cols": "functional"})

    # Prepare data for plotting functional clusters
    df_plot_melt = func_groups.melt("index", var_name='cols', value_name='vals').dropna().sort_values(by="index",
                                                                                                      ascending=True)

    # Adjust the index for plotting
    df_plot_melt["index"] = df_plot_melt["index"].apply(lambda x: x.split("_")[1])

    return df_plot_melt


def read_files_cdhit(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split("/")[-1]
        i = i.split(".cdhit")[0]
        #i = i.split("_")[0]
        dic[i] = element
    return dic

def read_files_interpro(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split("/")[-1]
        i = i.split(".csv")[0]
        #i = i.split(".tsv")[0]
        dic[i] = element
    return dic

def plot_all_sp(path_cdhit, path_interpro):
    dic_plot = {}
    for key_, path_ in read_files_interpro(path_interpro).items():
        for key, path in read_files_cdhit(path_cdhit).items():

            if key_ == key.split("_")[0]:
                dic_plot[key] = get_functional_cdhit_clusters(path, path_)
                dic_plot[key]["sp"] = key_

    df= pd.concat(dic_plot, axis=0).rename(columns={"cols":"categories", "sp":"species"})
    return df

df_plot2= plot_all_sp(path_cdhit, path_interpro).sort_values(by=["index", "vals"], ascending=[True, False])

""" Plot cdhit with  all species """

#colors={'#8c510a':'carpe','#bf812d': 'kbiala','#dfc27d': 'HIN','#f6e8c3' : 'trepo', '#80cdc1':'spiro','#35978f' : 'wb','#01665e': 'muris'}
sns.set_style("whitegrid", {'axes.grid': False})

colors = ['#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e']
pal= sns.set_palette(sns.color_palette(colors[::-1])) #reverse
#hue_order= ["carpe","kbiala","HIN","trepo", "spiro","wb", "muris"]

ax = sns.lineplot(y='vals', x="index", data=df_plot2, hue="species",palette=pal)
ax.set_xticklabels(['70%', '75%', '80%', '85%', '90%', '95%', '100%'])
ax.legend(title='species', loc='upper left',labels=['H. inflata', 'K. bialata',  'C. membranifera','Trepomonas pc1', 'S. salmonicida', 'G. intestinalis', 'G. muris'])

ax.set_xlabel("Sequence identity")
ax.set_ylabel("Functional clusters")
ax.axes.set_ylim(0,None)

plot = ax.get_figure()
#plot.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1b_all_sp.png', format="png", dpi=1200)
plot.show()
