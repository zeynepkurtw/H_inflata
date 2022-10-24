import pandas as pd
from Bio import SeqIO
import glob

#path = '/Users/zeyku390/PycharmProjects/H.inflata/resource/3_interproscan/*.tsv'
path = '/Users/zeyku390/PycharmProjects/H.inflata/resource/4_eggnog/eggnog/*.annotations'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[1]
    i = i.split("/")[-1]
    sp_dic[i] = element

