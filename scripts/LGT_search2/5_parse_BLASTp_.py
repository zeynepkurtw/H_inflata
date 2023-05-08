import pandas as pd
import glob
from ete3 import NCBITaxa


try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "output/4_BLASTp/concat/*.blastp"

    out_file = "data/LGT_search/blast_out/{}_filtered.blastp"

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
filtered = {}

for element in blast_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    filtered[i] = element

for key, values in filtered.items():

    df = pd.read_csv(values, sep='\t', header=None)
    #groupby query
    grouped = df.groupby(df.columns[0])

    #get top 100 unique sseqid
    result = grouped.apply(filter_unique_values)

    #make taxid str so that multiple taxids can be read
    result[15] = result[15].astype(str)
    result['kingdom'] = result[15].apply(get_kingdom_name)

    filtered[key] = result

for key, values in filtered.items():
    values.to_csv(out_file.format(key), sep="\t", index=False, header=None)





