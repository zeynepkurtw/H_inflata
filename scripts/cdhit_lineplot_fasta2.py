import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns

"""
input= Cdhit.fasta
output= df_func, df_hypo, df_id, df_plot, line_plot
"""

"read files from directory"
#path_cdhit = '/Users/zeyku390/PycharmProjects/H.inflata/output/2_cdhit/*.cdhit.clstr'
path_cdhit = '/Users/zeyku390/PycharmProjects/H.inflata/output/2_cdhit/cdhit_slow/*.cdhit'
path_interpro= "/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/HIN_interpro_annot.csv"

def read_files(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split(".fasta")[0]
        i = i.split("_")[-1]
        dic[i] = element
    return dic

""" get id lists from cdhit fasta
    return as dataframe"""
def get_id(fasta_path):
    id= []
    with open(fasta_path, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
    return pd.DataFrame(id)

"""" merge id lists with interpro annotation"
"""
dic_id={}
dic_ann={}
df_ann= pd.read_csv(path_interpro, sep="\t", header=None).drop_duplicates(0) #get unique genes
for key, values in read_files(path_cdhit).items():
    dic_id[key] = get_id(values)
    dic_ann[key]= pd.merge(dic_id[key], df_ann)


"Functional and hypothetical proteins"
def anti_join(x, y):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    # ans= ans.reset_index()
    return ans


dic_hypo = {}
dic_func = {}
for key, values in dic_ann.items():
    dic_func[key] = values
    dic_hypo[key] = anti_join(dic_id[key], dic_func[key])

"concatanate dic elements in to a dataframe"
df_id=pd.concat(dic_id, axis=1)
df_id.columns= df_id.columns.droplevel(1) #remove multindex columns

df_hypo=pd.concat(dic_hypo, axis=1)
df_hypo.columns= df_hypo.columns.droplevel(1) #remove multindex columns

df_func=pd.concat(dic_func, axis=1)
df_func.columns= df_func.columns.droplevel(1) #remove multindex columns

"get value counts from dataframe"
def get_counts(df):
    df_melt=df.melt( var_name='cols', value_name='vals').dropna().drop_duplicates()
    df_count=df_melt["cols"].value_counts().reset_index()
    return df_count


id_groups = get_counts(df_id).rename(columns={"cols": "id"})
hypo_groups = get_counts(df_hypo).rename(columns={"cols": "hypo"})
func_groups = get_counts(df_func).rename(columns={"cols": "func"})
#df_plot = pd.concat([id_groups, hypo_groups, func_groups], axis=0)
df_plot = pd.concat([id_groups, func_groups], axis=0)

df_plot_melt= df_plot.melt("index", var_name='cols', value_name='vals').dropna().sort_values(by="index",ascending=True)

"Plot cdhit results"

sns.set_style("whitegrid", {'axes.grid': False})
pal = sns.dark_palette("#69d", n_colors=2, reverse=True)

#ax = sns.lineplot(x='vals', y="index", data=df_plot_melt, hue="index", style="cols")
#ax = sns.lineplot(x='vals', y="index", data=df_plot_melt, hue="cols")
ax = sns.lineplot(y='vals', x="index", data=df_plot_melt, hue="cols", palette=pal)
ax.set_xticklabels(['70%', '75%', '80%', '85%', '90%', '95%', '100%'])
ax.set_xlabel("Sequence identity")
ax.set_ylabel("Number of clusters")
ax.legend(title='Clusters', loc='upper left', labels=['unique', 'functional'])

plot = ax.get_figure()
plot.show()