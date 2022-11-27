import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;

path_OG_count = "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
path_HOG= "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"


" Read files "
df_count = pd.read_csv(path_OG_count, sep="\t", header='infer', engine="python")
df = pd.read_csv(path_HOG, sep="\t", header='infer', engine="python")
df = pd.merge(df, df_count, left_on="OG", right_on="Orthogroup")

" 1:1:1:M single copy in parasites many copy in HIN "
df_1m = df[(df.HIN_y >= 1) & (df.iloc[:, [12,13,14,15,16,17]] == 1).all(1)]
plot1 = df_1m.iloc[:, [11,12,13,14,15,16,17]]
print(df_1m["HIN_y"].sum())



" Copy number distribution bar plot "
sns.set()
#sns.set(rc={'figure.figsize': (10, 5)})
sns.set_style("whitegrid", {'axes.grid': False})
#pal = sns.dark_palette("#69d", n_colors=1, reverse=True)
colors = ["#2A4A5C"]
pal = sns.set_palette(sns.color_palette(colors))

ax= sns.histplot(data=plot1['HIN_y'], discrete=True, palette=pal)
ax.set_xlabel('H. inflata',  style='italic')
ax.set_ylabel('Number of genes')

plt.xticks(range(1, 20))
#plt.title("m:1 H. inflata many copy other species single copy")
plt = ax.get_figure()
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1d.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()