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

    out_file = "data/LGT_search/blast_out/hin_trepo_cat_taxid.csv"
    out_file_xlsx = "data/LGT_search/blast_out/hin_trepo_cat_taxid.xlsx"

# read blast file top hits
columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
           "evalue", "bitscore", "qlen", "slen", "stitle", "staxids"]
df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
df= df.dropna(subset= "staxids", how ="any", axis=0)
df = df.drop_duplicates(subset="sseqid")

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

lineage = {}
names = {}
for taxon in df.loc[:, "staxids"]:
    try:

        if ';' in str(taxon):
            taxon = str(taxon).split(';')

            for part in taxon:
                lineage[part] = ncbi.get_lineage(part)
                names[part] = ncbi.get_taxid_translator(lineage[part])
        else:

            lineage[taxon] = ncbi.get_lineage(taxon)
            names[taxon] = ncbi.get_taxid_translator(lineage[taxon])

    except Exception as ex:
        print(ex)
        print(taxon)

taxid = pd.DataFrame(names).T.reset_index().rename(columns={'index': 'staxids'})
blast_out = pd.merge(df, taxid)

blast_out.to_csv(out_file, sep='\t', encoding='utf-8', index=False)
blast_out.to_excel(out_file_xlsx, index=False)