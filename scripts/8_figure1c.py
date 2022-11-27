import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import itertools

path_OG_count = "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
path_singletons = "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups_UnassignedGenes.tsv"

"OG with at least two genes"
df_count = pd.read_csv(path_OG_count, sep="\t", header='infer')
df_count = df_count.set_index("Orthogroup").sort_values(by="Total", ascending=False)
df_count.loc[df_count["Total"] > 1, "Type"] = "OG"

"Single genes which are excluded from the OG"
df_sing = pd.read_csv(path_singletons, sep="\t", header='infer')
df_sing = df_sing.set_index("Orthogroup").fillna(0)
df_sing["Total"] = 1
df_sing = df_sing.applymap(lambda x: 1 if isinstance(x, str) == True else x)

"Concatanate OG and singleton dataframes"
df_count_s = pd.concat([df_count, df_sing], axis=0)
df_count_s.loc[df_count_s["Total"] > 1, "Type"] = "OG"
df_count_s.loc[df_count_s["Total"] == 1, "Type"] = "singleton"

df_count_s = df_count_s.rename(
    columns={"HIN": "H. inflata", "spiro": "S. salmonicida", "trepo": "Trepomonas pc1",
             "wb": "G. intestinalis", "carpe": "C. membranifera", "kbiala": "K. bialata", "muris": "G. muris", })

# df_count_s=[["H. inflata","Trepomonas pc1" , "S. salmonicida", "G. intestinalis", "G. muris", "C. membranifera", "K. bialata"]]

sp_list = ["H. inflata", "G. muris", "S. salmonicida", "Trepomonas pc1",
           "G. intestinalis", "C. membranifera", "K. bialata"]

def count_types(df, sp):
    #df = df_count[[sp, "Type"]]
    # df= df[df[sp] != 0]
    other_cols = [i for i in sp_list if i != sp]
    df_ss = df[(df[sp] > 1) & (df[other_cols] == 0).all(axis=1)]
    df_ss["Type"] = "species-specific"

    df_sing= df[(df[sp] == 1) & (df[other_cols] == 0).all(axis=1)]

    df_new=pd.concat([df_ss, df_sing], axis=0)
    df_new = df_new["Type"].value_counts().reset_index()
    df_new["sp"] = sp
    return df_new


dic_count = {}
for sp in sp_list:
    dic_count[sp] = count_types(df_count_s, sp)
df_plot = pd.concat(dic_count, axis=0).sort_values(by=["index", "Type"], ascending=[True, False])

""" Bar plot    """

sns.set_style("whitegrid", {'axes.grid': False})
#pal = sns.dark_palette("#69d", n_colors=2, reverse=True)
colors = ["#7A93A2", "#2A4A5C"]

pal = sns.set_palette(sns.color_palette(colors))

order = ['C. membranifera',  'K. bialata', 'H. inflata', 'Trepomonas pc1', 'S. salmonicida', 'G. intestinalis', 'G. muris']
ax = sns.barplot(x='sp', y='Type', hue='index', data=df_plot, palette=pal, order=order )

hatches = [ "", "." ]
# Loop over the bars
for bars, hatch in zip(ax.containers, hatches):
    # Set a different hatch for each group of bars
    for bar in bars:
        bar.set_hatch(hatch)


ax.legend().set_title("Group categories")
ax.set_xlabel('')
ax.set_ylabel('Number of categories')
plt.xticks(rotation=45, ha="right" )
# ax.set_xticklabels( labels="index", rotation=45, style='italic',ha="right")

plot = ax.get_figure()
plot.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1c.png', format="png", bbox_inches='tight',
             dpi=1200)
plot.show()
