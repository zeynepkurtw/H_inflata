import pandas as pd
from Bio import SeqIO
import glob

"""
Add interproscan entries
input: sp_interpro.tsv
output: sp_interpro_annot.csv

"""


path = '/Users/zeyku390/PycharmProjects/H_inflata/resource/3_interproscan/*.tsv'
#path = '/Users/zeyku390/PycharmProjects/H_inflata/resource/4_eggnog/eggnog/*.annotations'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element


dic_annot= {}
for key, value in sp_dic.items():
    dic_annot[key] = pd.read_csv(value, sep="\t", header=None, names=list(range(0, 15)), engine='python', quoting=3)[[0,11,12]]
    dic_annot[key].dropna().drop_duplicates().to_csv(f"/Users/zeyku390/PycharmProjects/H_inflata/jupyter/data/interpro_annot/{key}.csv",
                                 sep="\t", header=None, index=False)


