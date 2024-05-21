import pandas as pd
from upsetplot import UpSet
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

path_OG_count = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
out_1 = 'plots/Orthofinder_7sp/upset1_comb.png'


"OG with at least two genes"
df_count = pd.read_csv(path_OG_count, sep="\t", header='infer')
df_count = df_count.set_index("Orthogroup").sort_values(by="Total", ascending=False)
df_count.loc[df_count["Total"] > 1, "Type"] = "OG"


df_count = df_count.rename(columns={ "HIN": "H. inflata",
                                        "trepo": "Trepomonas pc1",
                                        "spiro": "S. salmonicida",
                                        "wb": "G. intestinalis",
                                        "muris": "G. muris",
                                        "kbiala": "K. bialata",
                                        "carpe": "C. membranifera",

                                         })
df_count =df_count[["C. membranifera", "K. bialata","H. inflata", "Trepomonas pc1", "S. salmonicida",
                        "G. intestinalis", "G. muris", "Total", "Type"]]



"Upset plot for both OG and singletons marked different colors"
df_stack = df_count.set_index(df_count["C. membranifera"] >= 1). \
    set_index(df_count["K. bialata"] >= 1, append=True). \
    set_index(df_count["H. inflata"] >= 1, append=True). \
    set_index(df_count["Trepomonas pc1"] >= 1, append=True). \
    set_index(df_count["S. salmonicida"] >= 1, append=True). \
    set_index(df_count["G. intestinalis"] >= 1, append=True). \
    set_index(df_count["G. muris"] >= 1, append=True)


""" plot upset0 """
upset0 = UpSet(df_stack,
               min_subset_size=20,
               min_degree=2,
               show_counts=True,
               sort_by="cardinality",
               facecolor="darkgray",
               intersection_plot_elements=5)

upset0.style_subsets(absent=["K. bialata", "C. membranifera"],
                         min_degree=5,
                         facecolor="salmon",
                         label="Diplomonads")

upset0.style_subsets(absent=["K. bialata", "C. membranifera", "S. salmonicida", "G. intestinalis", "G. muris"],
                     facecolor="darkseagreen",
                     label="2ndary Free-livings")

upset0.style_subsets(absent=["S. salmonicida", "G. intestinalis", "G. muris"],
                     min_degree=4,
                     facecolor="steelblue",
                     label="Free-livings")

pal={"C. membranifera": "steelblue",
                            "H. inflata": "darkseagreen",
                            "K. bialata": "steelblue",
                            "Trepomonas pc1": "darkseagreen",
                            "S. salmonicida": "salmon",
                            "G. intestinalis": "salmon",
                            "G. muris": "salmon",
     "Total": "lightgray"
                            }

upset0.add_catplot(value="C. membranifera",
                   kind="box",
                   color="steelblue",
                     )

upset0.add_catplot(value="K. bialata",
                   kind="box",
                   color="steelblue",
                     )

upset0.add_catplot(value="Trepomonas pc1",
                   kind="box",
                  color="darkseagreen",
                     )

upset0.add_catplot(value="H. inflata",
                   kind="box",
                  color="darkseagreen",
                     )


upset0.add_catplot(value="S. salmonicida",
                   kind="box",
                   color=   "salmon",
                     )

upset0.add_catplot(value="G. intestinalis",
                   kind="box",
                   color=   "salmon",
                     )

upset0.add_catplot(value="G. muris",
                   kind="box",
                  color=  "salmon",
                     )


plt.title("Count by type")
sns.set_style("whitegrid", {'axes.grid': False})
upset0.plot()

plt.savefig(out_1, format="png", bbox_inches='tight', dpi=800)
plt.show()
