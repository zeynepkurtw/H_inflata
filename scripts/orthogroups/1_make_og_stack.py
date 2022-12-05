import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt"

    out_file =  "data/orthogroups/og_stack.csv"


df = pd.read_csv(og_file, header=None, dtype=str, delim_whitespace=True)
df = df.apply(lambda x: x.str.replace(":", ""))
df = df.rename(columns={0: "OG"}).set_index("OG")
df = df.stack().reset_index().drop(columns=["level_1"])
df = df.rename(columns={0:"id"})

df.to_csv(out_file, sep="\t", index=False)
