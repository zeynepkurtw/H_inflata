import pandas as pd
import glob
try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_file = "data/orthogroups/og_ann.csv"
    ipr_files= "data/superfamily/signature_iprs_*"

    out_file =  "data/superfamily/og_ann_ipr_LCA.csv" #lrr, cystine, ankyrin

family = {}
list_files = glob.glob(ipr_files)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    j = i.split("_")[-1]

    family[j] = element

#get Leucine, cysteine and ankyrin families by signature IPR

for key, value in family.items():
    og_ann = pd.read_csv(og_ann_file,  header="infer", sep="\t")
    ipr = pd.read_csv(value,  header="infer", sep="\t")
    ipr["family"] = key

    family[key] = pd.merge(og_ann, ipr, on="ipr")

    print(f"{key} id ",len(family[key]["id"].drop_duplicates()) )
    print(f"{key} ipr ",len(family[key]["ipr"].drop_duplicates()) )
    print(f"{key} id-ipr",len(family[key].groupby(["id", "ipr"]).size() ) )
    print("")

#og_ann_ipr_LCA = pd.DataFrame.from_dict(family)
og_ann_ipr_LCA = pd.concat(family, axis=0).sort_values("id")
og_ann_ipr_LCA.to_csv(out_file.format(key), index=False, sep="\t")


