import pandas as pd

"""
Signature IPRs

NOTE= Curate signature IPR tsv file

"""

import pandas as pd
import glob

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    ipr_files = "resource/6_superfamily/*.tsv"

    out_file =  "data/superfamily/signature_iprs_{}.tsv" #lrr, cystine, ankyrin

family = {}
list_files = glob.glob(ipr_files)
for element in list_files:
    try:
        i = element.split(".")[0]
        i = i.split("/")[-1]
        j = i.split("-")[0]
        k = i.split("-")[1]
        family[i, j, k] = element
    except Exception as ex:
        print(ex)


for key, value in family.items():

    #read ipr file
    df=pd.read_csv(value, header="infer", sep="\t")

    #get IPRs
    ipr=df[df["Accession"].str.contains("IPR", na=False)].sort_values(by="Name")

    # get family IPRs
    ipr_1=ipr[ipr["Name"].str.contains(key[1], na=False , case= False)]
    ipr_2=ipr[ipr["Name"].str.contains(key[2], na=False , case= False)]

    all_iprs = pd.concat([ipr_1, ipr_2], axis=0).drop_duplicates(subset=["Accession"])[["Accession"]]
    all_iprs = all_iprs.rename(columns= {"Accession": "ipr"})

    #Filter IPRs in Description & Name
    print("signature_iprs_"+ key[0],len(all_iprs))

    all_iprs.to_csv(out_file.format(key[0]), index=False, sep="\t")