import pandas as pd
from upsetplot import UpSet
from upsetplot import generate_counts, plot
from matplotlib import pyplot as plt
import seaborn as sns

path_OG_count = "output/1_orthofinder/Results_Nov28/Orthogroups/Orthogroups.GeneCount.tsv"
path_singletons = "output/1_orthofinder/Results_Nov28/Orthogroups/Orthogroups_UnassignedGenes.tsv"
out_0 = 'plots/Orthofinder_10sp/upset0_10sp.png'
out_1 = 'plots/Orthofinder_10sp/upset1_10sp.svg'

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

df_count_s = df_count_s.rename(columns={"gruberi": "N. gruberi",
                                        "exilis": "M. exilis",
                                        "trico": "T. vaginalis",
                                        "carpe": "C. membranifera",
                                        "kbiala": "K. bialata",
                                        "HIN": "H. inflata",
                                        "trepo": "Trepomonas pc1",
                                        "spiro": "S. salmonicida",
                                        "wb": "G. intestinalis",
                                        "muris": "G. muris"})
# df_count_s =df_count_s[["H. inflata", "K. bialata" ]]

"Upset plot for both OG and singletons marked different colors"
df_stack = df_count_s.set_index(df_count_s["N. gruberi"] >= 1). \
    set_index(df_count_s["M. exilis"] >= 1, append=True). \
    set_index(df_count_s["T. vaginalis"] >= 1, append=True). \
    set_index(df_count_s["C. membranifera"] >= 1, append=True). \
    set_index(df_count_s["K. bialata"] >= 1, append=True). \
    set_index(df_count_s["H. inflata"] >= 1, append=True). \
    set_index(df_count_s["Trepomonas pc1"] >= 1, append=True). \
    set_index(df_count_s["S. salmonicida"] >= 1, append=True). \
    set_index(df_count_s["G. intestinalis"] >= 1, append=True). \
    set_index(df_count_s["G. muris"] >= 1, append=True)

""" plot upset0 """
upset0 = UpSet(df_stack.sort_values(by="Total", ascending=True),
               # max_degree= 1,
               min_subset_size=10,
               intersection_plot_elements=0,
               show_counts=True,
               # sort_categories_by="cardinality",
               # sort_by="cardinality")  # disable the default bar chart
               )

pal = sns.dark_palette("#69d", n_colors=2, reverse=False)
upset0.add_stacked_bars(by="Type",
                        colors=pal,
                        title="Count by type",
                        elements=5)

upset0.style_subsets(max_degree=1,
                     facecolor="white",
                     edgecolor="black",
                     label="Species-specific")

sns.set_style("whitegrid", {'axes.grid': False})
upset0.plot()
plt.suptitle("OG Upset plot")
#plt.savefig(out_0, format="png", bbox_inches='tight', dpi=1200)
#plt.show()

""" plot upset1: Filtered"""

upset1 = UpSet(df_stack,
               min_degree=2,
               min_subset_size=10,
               show_counts=True,
               sort_categories_by=None,
               intersection_plot_elements=10,
               sort_by="degree",
               # facecolor="#2A4A5C",
               facecolor="#3A5E7A")  # disable the default bar chart

upset1.style_subsets(present=["H. inflata"],
                    min_degree=2,
                    facecolor="#a50026",
                     label="OGs shared by Genes shared by H. inflata")


sns.set_style("whitegrid", {'axes.grid': False})
upset1.plot()

# plt.suptitle(" Orthogroups Upset plot ")
plt.savefig(out_1, format="svg", bbox_inches='tight', dpi=1200)
plt.show()
