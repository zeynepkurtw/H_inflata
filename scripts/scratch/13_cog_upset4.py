import pandas as pd
from upsetplot import UpSet
from upsetplot import from_memberships

from matplotlib import pyplot as plt
import seaborn as sns
from utils import get_og_stack, get_og_count_sing

"""
OG gene list stack
"""

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann = "data/orthogroups/og_ann.csv"

    out_file = "plots/upset_cog4.png"


df = pd.read_csv(og_ann, sep="\t", header='infer',
                 dtype={"carpe": int, "kbiala": int, "HIN": int, "trepo": int, "wb": int, "spiro": int, "muris": int,
                        "Total": int, "OG": str, "id": str, "ann_f": str, "ann_inter": str, "ann_egg": str, "ipr": str,
                        "COG": str, "KEGG_KOs": str, "LGT": str, "Type": str})
df = df[["OG","id",
         "carpe",
         "kbiala",
         "HIN",
         "trepo",
         "spiro",
         "wb",
         "muris",
         "COG"]]

df = df.dropna(subset=['COG'])
df = df.drop_duplicates()
df["COG"] = df["COG"].str.split(",")
df = df.explode("COG").reset_index(drop=True)

genes_by_cog = from_memberships(df.COG.str.split(','), data=df)


UpSet(genes_by_cog, min_subset_size=15, show_counts=True).plot()

plt.savefig(out_file, format="png", dpi=1200, bbox_inches='tight')
plt.show()