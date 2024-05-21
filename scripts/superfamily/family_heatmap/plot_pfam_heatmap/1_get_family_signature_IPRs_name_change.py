import pandas as pd
import glob

"""
Extract IPR entries from SearchResults that are downloaded from interproscan website
"""

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    ipr_files = "resource/6_superfamily/name_change/SearchResults*"

    out_file = "data/superfamily/1_signature_iprs/plot_pfam_heatmap/signature_iprs_{}.tsv"  # top 7 sueprfamily

family = {}
list_files = glob.glob(ipr_files)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    j = i.split("-")[1]
    family[i, j] = element

ipr_dict ={}
for key, value in family.items():
    # read search-result file
    df = pd.read_csv(value, header="infer", sep="\t")
    # get IPRs
    ipr = df[df["Accession"].str.contains("IPR", na=False)].rename(columns={"Accession": "ipr"})
   # ipr_dict[key[1]] = ipr[ipr["Name"].str.contains(key[1], na=False, case=False)]
    ipr_dict[key[1]] = ipr.drop_duplicates(subset=["ipr"])[["ipr"]]

    print("j = ", key[1])
    print("# ipr = ", len(ipr_dict[key[1]]))

    ipr_dict[key[1]].to_csv(out_file.format(key[1]), index=False, sep="\t")
