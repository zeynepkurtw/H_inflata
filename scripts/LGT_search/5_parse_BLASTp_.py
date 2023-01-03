import pandas as pd
import glob
from Bio import SeqIO
from Bio import Entrez

try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "output/4_BLASTp/hin_trepo_cat.blastp"
    taxa_file = "data/LGT_search/blast_out/hin_trepo_cat_taxa_.blastp"

    out_file = "data/LGT_search/blast_out/hin_trepo_cat_taxa_bacteria.blastp"


df = pd.read_csv(blast_file, sep='\t', header=None)
df_thit = df.drop_duplicates(subset=0, keep="first")



"""hit_list = []
for i in df.loc[:, "qseqid"]:
    if i == "Eukaryota"
        hit_list.append(hit)
"""


