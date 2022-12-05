import pandas as pd
from Bio import SeqIO
import glob


path = 'resource/4_eggnog/eggnog/*.annotations'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    i = i.split("_")[0]
    sp_dic[i] = element


dic_ann= {}
for key, value in sp_dic.items():
    dic_ann[key] = pd.read_csv(value, sep="\t",  header="infer", skiprows=3)[["#query_name","COG cat", "KEGG_KOs", "eggNOG annot"]]
    dic_ann[key] = dic_ann[key].drop_duplicates().rename(columns={"#query_name": "id", "COG cat":"COG", "eggNOG annot":"ann_egg" })
    dic_ann[key].to_csv(f"data/eggnog_ann/{key}.csv", sep="\t", index=False)

pd.concat(dic_ann, axis=0).dropna().drop_duplicates().to_csv("data/eggnog_ann/ann_egg_cat.csv", sep="\t", index=False)
