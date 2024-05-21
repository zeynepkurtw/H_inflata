import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns

"""
input= Cdhit.fasta
output= df_func, df_hypo, df_id, df_plot, line_plot
"""

"read files from directory"
#path_cdhit = '/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/*.cdhit.clstr'
path_cdhit = '/output/2_cdhit/cdhit_fast/*.cdhit.clstr'
path_interpro= "/jupyter/data/HIN_interpro_annot.csv"

def read_files(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split(".fasta")[0]
        i = i.split("_")[-1]
        dic[i] = element
    return dic



"""" merge id lists with interpro annotation"
"""
dic_id={}
dic_ann={}
dic_clstr={}
df_ann= pd.read_csv(path_interpro, sep="\t", header=None)
for key, values in read_files(path_cdhit).items():
    dic_clstr[key]= pd.read_csv(values, sep="\t", header=None )
    #dic_id[key] = get_id(values)
    dic_ann[key]= pd.merge(dic_id[key], df_ann)


