import pandas as pd
from ete3 import NCBITaxa
import warnings
warnings.filterwarnings('ignore')

try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "output/4_BLASTp/hin_trepo_cat.blastp"

    out_file = "data/LGT_search/blast_out/hin_trepo_cat_stats.csv"
    out_file_xlsx = "data/LGT_search/blast_out/hin_trepo_cat_stats.xlsx"

# read blast file top hits
columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
           "evalue", "bitscore", "qlen", "slen", "stitle", "staxids"]
df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
df= df.dropna(subset= "staxids", how ="any", axis=0)[0:500]



lineage = {}
names = {}
ncbi = NCBITaxa()
for taxon in df.loc[:, "staxids"]:
    if ';' in str(taxon):
        taxon = str(taxon).split(';')

        for part in taxon:
            lineage[part] = ncbi.get_lineage(part)
            names[part] = ncbi.get_taxid_translator(lineage[part])

    else:
        lineage[taxon] = ncbi.get_lineage(taxon)
        names[taxon] = ncbi.get_taxid_translator(lineage[taxon])


taxid = pd.DataFrame(names).T.reset_index().rename(columns={'index': 'staxids'})
blast_out = pd.merge(df, taxid)

gb = blast_out.groupby(["qseqid", "sseqid"])

for gene_id, group in gb:
    # get the top 100 hits
    top100 = group.head(100)
    # calculate the percentage of different taxons
    stats = top100['staxids'].value_counts() / len(top100)
    # add the stats as a new column
    data.loc[data['gene_id'] == gene_id, 'stats'] = stats





blast_out.to_csv(out_file, sep='\t', encoding='utf-8', index=False)
blast_out.to_excel(out_file_xlsx, index=False)
