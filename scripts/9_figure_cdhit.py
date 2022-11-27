import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns
import matplotlib as plt

"""
Cdhit line plot
input= cdhit.fasta
output= df_func, df_hypo, df_id, df_plot, line_plot
"""

"read files from directory"
#path_cdhit = '/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/*.cdhit.clstr'
path_cdhit = "/Users/zeyku390/PycharmProjects/H_inflata/output/2_cdhit/cdhit_fasta/*.cdhit"
path_interpro= "/Users/zeyku390/PycharmProjects/H_inflata/jupyter/data/interpro_annot/*.csv"

def read_files_cdhit(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split(".cdhit")[0]
        i = i.split("/")[-1]
        dic[i] = element
    return dic

def read_files_interpro(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split(".")[0]
        i = i.split("/")[-1]
        dic[i] = element
    return dic

""" get id lists from cdhit fasta
    return as dataframe """
def get_id(fasta_path):
    id= []
    with open(fasta_path, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
    return pd.DataFrame(id)

dic_id={}
dic_ann={}
for i,j in read_files_cdhit(path_cdhit).items():
    for m,n in read_files_interpro(path_interpro).items():
        dic_id[i] = get_id(j)
        dic_ann[m] = pd.read_csv(n, sep="\t", header=None).drop_duplicates(0) #get unique genes

"""Functional and hypothetical proteins"""
def anti_join(x, y):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    # ans= ans.reset_index()
    return ans

dic_hypo = {}
dic_func = {}
for key, values in dic_ann.items():
    dic_func[key] = values
    dic_hypo[key] = anti_join(dic_id[key], dic_func[key])
