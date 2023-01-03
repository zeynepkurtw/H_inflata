import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa

try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "output/4_BLASTp/hin_trepo_cat.blastp"

    out_file = "data/LGT_search/blast_out/hin_trepo_cat_taxa_.blastp"

# read blast file
columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
           "evalue", "bitscore", "qlen", "slen", "stitle", "staxids"]
df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
#df = df.drop_duplicates(subset="sseqid")[1:10]


def check_LGT(taxon, LGT):
    ncbi = NCBITaxa()
    #ncbi.update_taxonomy_database()
    lineage = ncbi.get_lineage(taxon)
    names = ncbi.get_taxid_translator(lineage)

    if 'Eukaryota' not in [names[taxid] for taxid in lineage]:
        LGT = True
    return LGT


#get putative LGT list
protein_results = df.groupby('qseqid')

for query, df in protein_results:
    LGT = False
    taxa = df['staxids'].dropna().values.tolist()

    for taxon in taxa:
        if ';' in str(taxon):
            taxon = str(taxon).split(';')

            for part in taxon:
                LGT = check_LGT(part, LGT)

        else:
            LGT = check_LGT(taxon, LGT)

    if LGT == True:
        pd.DataFrame(taxa).to_csv(out_file, sep="\t")



"""protein_results = df.groupby('qseqid')

        for query, df in protein_results:
            eukaryote = True
            taxa = df['staxids'].dropna().values.tolist()

            for taxon in taxa:
                if ';' in str(taxon):
                    taxon = str(taxon).split(';')

                    for part in taxon:
                        eukaryote = True
                else:
                    eukaryote = False

        if eukaryote == False:
            to_csv()
"""
