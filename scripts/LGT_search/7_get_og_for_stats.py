import pandas as pd
import glob


try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    stats_file = "data/LGT_search/blast_out/*_filtered_stats.blastp"
    og_ann= "data/orthogroups/og_ann.csv"

    out_file = "data/LGT_search/blast_out/{}_og.blastp"



stats_files = glob.glob(stats_file)
stats_og = {}

for element in stats_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    stats_og[i] = element

for key, values in stats_og.items():
    df_stats = pd.read_csv(values, sep='\t', header="infer")
    df_og_ann= pd.read_csv(og_ann, sep="\t", header="infer")[["OG", "id"]]
    stats_og[key] = pd.merge(df_stats, df_og_ann)[["OG", "id", "taxa", "perc_taxa"]].sort_values(by=["OG", "id"])
    stats_og[key].to_csv(out_file.format(key), sep="\t", index=False)