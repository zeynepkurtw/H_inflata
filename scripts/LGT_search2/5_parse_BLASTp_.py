import pandas as pd
import glob
from ete3 import NCBITaxa


try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "output/4_BLASTp/concat/*.blastp"
    #taxa_file = "data/LGT_search/blast_out/hin_trepo_cat_taxa_.blastp"

    out_file = "data/LGT_search/blast_out/{}_stats.blastp"

def filter_unique_values(group):
    unique_values = group[group.columns[1]].unique()[:100]
    return group[group[group.columns[1]].isin(unique_values)]

ncbi = NCBITaxa()


def get_kingdom_name(taxon_id):
    if ';' in taxon_id:
        taxon_ids = taxon_id.split(';')
        taxon_id = taxon_ids[0]
    taxon_id = pd.to_numeric(taxon_id, errors='coerce')
    if pd.isna(taxon_id):
        return None

    try:
        lineage = ncbi.get_lineage(int(taxon_id))
        if lineage is not None:
            ranks = ncbi.get_rank(lineage)
            for key , values in ranks.items():
                if values == "superkingdom":
                    names = ncbi.get_taxid_translator([key])[key]
                    return names

    except:
        return None
    return None



blast_files = glob.glob(blast_file)
for file in blast_files:
    df = pd.read_csv(file, sep='\t', header=None)
    # get unique sseqid (already sorted by bitscore)
    #df = df.drop_duplicates(subset=1, keep="first")

    #groupby query
    grouped = df.groupby(df.columns[0])
    #get top 100 unique sseqid
    result = grouped.apply(filter_unique_values)
    #make taxid str so that multiple taxids can be read
    result[15] = result[15].astype(str)
    result['kingdom'] = result[15].apply(get_kingdom_name)
