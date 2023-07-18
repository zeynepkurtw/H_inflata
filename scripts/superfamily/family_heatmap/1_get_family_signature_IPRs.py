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
    ipr_files = "resource/6_superfamily/SearchResults*"

    out_file = "data/superfamily/1_signature_iprs/signature_iprs_{}.tsv"  # lrr, cystine, ankyrin

family = {}
list_files = glob.glob(ipr_files)
for element in list_files:
    try:
        i = element.split(".")[0]
        i = i.split("/")[-1]
        j = i.split("-")[1]
        k = j.split(" ")[0]
        family[i, j, k] = element

    except Exception as ex:
        i = element.split(".")[0]
        i = i.split("/")[-1]
        j = i.split("-")[1]
        k = j.split(" ")[0]
        l = k.split(" ")[0]
        family[i, j, k, l] = element

ipr_dict ={}
for key, value in family.items():
    # read ipr file
    df = pd.read_csv(value, header="infer", sep="\t")
    # get IPRs
    ipr = df[df["Accession"].str.contains("IPR", na=False)].rename(columns={"Accession": "ipr"})
    ipr_dict[key[2]] = ipr[ipr["Name"].str.contains(key[2], na=False, case=False)]
    ipr_dict[key[2]] = ipr_dict[key[2]].drop_duplicates(subset=["ipr"])[["ipr"]]

    print("j = ", key[1])
    print("k = ", key[2])
    print("# ipr = ", len(ipr_dict[key[2]]))

    ipr_dict[key[2]].to_csv(out_file.format(key[1]), index=False, sep="\t")
