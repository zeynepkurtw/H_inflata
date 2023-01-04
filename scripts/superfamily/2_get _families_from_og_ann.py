import pandas as pd
import glob
try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_file = "data/orthogroups/og_ann.csv"
    ipr_files= "data/superfamily/signature_iprs_*"

    out_file =  "data/superfamily/og_ann_family_{}.csv" #lrr, cystine, ankyrin

family = {}
list_files = glob.glob(ipr_files)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    j = i.split("_")[-1]

    family[j] = element

#get families by signature IPR

for key, value in family.items():
    og_ann = pd.read_csv(og_ann_file,  header="infer", sep="\t")
    ipr = pd.read_csv(value,  header="infer", sep="\t")

    df= pd.merge(og_ann, ipr, on="ipr")
    df["family"] = key
    df.to_csv(out_file.format(key), index=False, sep="\t")

#get families by annotation

