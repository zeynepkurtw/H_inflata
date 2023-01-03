# extract taxid, pi, stitle
# top hit
# divide hits into prokaryotes and eukaryotes
# detect putative LGT

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

    out_file = "data/LGT_search/blast_out/hin_trepo_cat_taxa2.blastp"

# read blast file
columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
           "evalue", "bitscore", "qlen", "slen", "stitle", "staxids"]
df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)

#take only top hits
df = df.drop_duplicates(subset="qseqid")#[1:100]

# make list from hits
hit_list = []
for hit in df.loc[:, "sseqid"]:
    hit_list.append(hit)

# make phylum df from entrez
Entrez.email = "zeynep.kurt@icm.uu.se"
taxa = []
species_name = []

for i in hit_list:

    handle = Entrez.efetch(db="protein", id=i, retmode="xml")
    record = Entrez.read(handle)
    taxa.append(record[0]["GBSeq_taxonomy"])
    species_name.append(record[0]["GBSeq_organism"])


phylum_df = pd.DataFrame(list(zip(hit_list, species_name)),
                         columns=["sseqid", "species"])

# merge blast with entrez info
blast_out = pd.merge(df, phylum_df)
blast_out.to_csv(out_file, sep='\t', encoding='utf-8', index=False)
