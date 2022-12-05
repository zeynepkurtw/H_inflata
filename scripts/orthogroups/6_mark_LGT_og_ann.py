import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann = "data/orthogroups/og_ann.csv"
    trepo_file = "data/LGT/trepo_lgt.csv"

    out_file = "data/orthogroups/og_ann_lgt.csv"


og_ann_ = pd.read_csv(og_ann, sep="\t", header="infer")
trepo = pd.read_csv(trepo_file, sep="\t", header=None)[[0, 1]].rename(columns={0: "OG", 1: "id"})
trepo["LGT"] = "trepo"

df = pd.merge(og_ann_, trepo, how="left").sort_values(by=["LGT", "OG", "Total"], ascending=[False, False, False])
df.to_csv(out_file, sep="\t", index=False)
