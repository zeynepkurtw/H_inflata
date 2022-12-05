from utils import *


try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_count_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
    trepo_file = "data/LGT/trepo_lgt.csv"

    outfile_csv = "data/og_groups/og_trepo_lgt.csv"
    outfile_xlsx = "data/og_groups/og_trepo_lgt.xlsx"


df = get_ann(get_og_genes(get_og_count(og_count_file)))
df = merge_trepo_lgt(trepo_file)

df.to_csv(outfile_csv, sep="\t", index=False)
df.to_excel(outfile_xlsx, index=False)
