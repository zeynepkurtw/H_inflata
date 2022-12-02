import pandas as pd
from upsetplot import UpSet
from upsetplot import generate_counts, plot
from matplotlib import pyplot as plt
import seaborn as sns

path_OG_count = "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
path_singletons = "/Users/zeyku390/PycharmProjects/H_inflata/output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups_UnassignedGenes.tsv"

"OG with at least two genes"
df_count = pd.read_csv(path_OG_count , sep="\t", header='infer')
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

df_count_s = df_count_s.rename(columns = {"carpe":"C. membranifera",
                                        "kbiala":"K. bialata",
                                          "HIN":"H. inflata",
                                        "trepo": "Trepomonas pc1",
                                          "spiro":"S. salmonicida",
                                          "wb":"G. intestinalis",
                                            "muris":"G. muris"})
#df_count_s =df_count_s[["H. inflata", "K. bialata" ]]

"Upset plot for both OG and singletons marked different colors"
df_stack = df_count_s.set_index(df_count_s["C. membranifera"] >= 1).\
    set_index(df_count_s["K. bialata"] >= 1, append=True).\
    set_index(df_count_s["H. inflata"] >= 1, append=True).\
    set_index(df_count_s["Trepomonas pc1"]>= 1, append=True).\
    set_index(df_count_s["S. salmonicida"] >= 1, append=True).\
    set_index(df_count_s["G. intestinalis"] >= 1, append=True).\
    set_index(df_count_s["G. muris"] >= 1, append=True)

""" plot upset0 """
upset0 = UpSet(df_stack.sort_values(by="Total", ascending=True),
              #max_degree= 1,
              min_subset_size=10,
              intersection_plot_elements=0,
              show_counts=True,
              #sort_categories_by="cardinality",
              #sort_by="cardinality")  # disable the default bar chart
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
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up0.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()




""" plot upset1: Filtered"""

upset1 = UpSet(df_stack,
              min_degree= 2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by=None,
             intersection_plot_elements=10,
              sort_by="degree",
               #facecolor="#2A4A5C",
               facecolor="#3A5E7A")  # disable the default bar chart


sns.set_style("whitegrid", {'axes.grid': False})
upset1.plot()

#plt.suptitle(" Orthogroups Upset plot ")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up1.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()

""" plot upset2: with singletons """

upset2 = UpSet(df_stack.sort_values(by="Total", ascending=True),
              max_degree= 1,
              #min_subset_size=10,
              intersection_plot_elements=0,
              show_counts=True,
              #sort_categories_by="cardinality",
              #sort_by="cardinality")  # disable the default bar chart
               )

pal = sns.dark_palette("#69d", n_colors=2, reverse=False)
upset2.add_stacked_bars(by="Type",
                       colors=pal,
                       title="Count by type",
                       elements=5)

upset2.style_subsets(max_degree=1,
                    #facecolor="white",
                    #edgecolor="black",
                    label="Species-specific")

sns.set_style("whitegrid", {'axes.grid': False})
upset2.plot()
#plt.suptitle("Singleton and species-specific OG Upset plot")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up2.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()

""" plot upset3: diplomonads marked """

upset3 = UpSet(df_stack,
              intersection_plot_elements=10,
              min_degree=2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by="cardinality",
              sort_by="cardinality")

upset3.style_subsets(absent=["K. bialata", "C. membranifera"],
                    min_degree=5,
                    facecolor="blue",
                    label="OGs shared by all diplomonads")

#{"HIN":"H. inflata", "muris":"G. muris", "S. salmonicida":"S. salmonicida", "trepo": "Trepomonas pc1", "G. intestinalis":"G. intestinalis", "carpe":"C. membranifera","kbiala":"K. bialata"})


upset3.plot()
#plt.suptitle("OGs shared by all diplomonads")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up3.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()


""" plot upset4: HIN and trepo marked """
upset4 = UpSet(df_stack,
              intersection_plot_elements=10,
              min_degree=2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by="cardinality",
              sort_by="cardinality")

upset4.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                    facecolor="red",
                    label="OGs shared by H. inflata and Trepomonas pc1")

upset4.plot()
#plt.suptitle("OGs shared by HIN and trepo")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up4.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()


""" plot upset5: FL marked """
upset5 = UpSet(df_stack,
              intersection_plot_elements=10,
              min_degree=2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by="cardinality",
              sort_by="cardinality")


upset5.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                    min_degree=4,
                    facecolor="green",
                    label="OGs shared by Free-living species")

upset5.plot()
#plt.suptitle("OGs shared by Free-living species")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up5.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()


""" plot upset6: HIN, trepo, carpe marked """
upset6 = UpSet(df_stack,
              intersection_plot_elements=10,
              min_degree=2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by="cardinality",
              sort_by="cardinality")


upset6.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris", "K. bialata"],
                    min_degree=3,
                    facecolor="green",
                    edgecolor="red",
                    label="OGs shared by by H. inflata, Trepomonas, C. membranifera")

upset6.plot()
#plt.suptitle("OGs shared by by H. inflata, Trepomonas, C. membranifera")
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up6.png', format="png", bbox_inches='tight', dpi=1200)
plt.show()


""" plot upset7: HIN, trepo, kbiala marked """
upset7 = UpSet(df_stack,
              intersection_plot_elements=10,
              min_degree=2,
              min_subset_size=10,
              show_counts=True,
              sort_categories_by="cardinality",
              sort_by="cardinality")

upset7.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris", "C. membranifera"],
                    min_degree=3,
                    facecolor="green",
                    edgecolor="red",
                     label="OGs shared by Genes shared by H. inflata, Trepomonas, K. bialata")

upset7.plot()
#plt.suptitle("OGs shared by Genes shared by H. inflata, Trepomonas, K. bialata")
#plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up7.svg', format="svg", bbox_inches='tight', dpi=1200)
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1e_up7.png', format="png", bbox_inches='tight', dpi=1200)

plt.show()