import pandas as pd
import glob
try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_file = "data/orthogroups/og_ann.csv"
    aa_perc_files = "data/superfamily/aa_percent-*"

    out_file =  "data/superfamily/og_ann_aa.csv" #lrr, cystine, ankyrin

family = {}
list_files = glob.glob(aa_perc_files)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    j = i.split("-")[-1]
    family[j] = element

#Add aa percent to each genes
dic_aa= {}
for key, value in family.items():
    og_ann = pd.read_csv(og_ann_file,  header="infer", sep="\t")
    aa = pd.read_csv(value,  header="infer", sep="\t")
    dic_aa[key]= pd.merge(og_ann, aa, on="id")


og_ann_aa = pd.concat(dic_aa, axis=0).dropna().drop_duplicates()
og_ann_aa.to_csv(out_file, sep="\t", index=False)
