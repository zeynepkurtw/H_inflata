import pandas as pd
from Bio import SeqIO
import glob

path = '/Users/zeyku390/PycharmProjects/data/proteome/*.fa'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element


def get_annot( fasta ):
    """
    Get a list of genes with annotations from fasta.
    :param old: FASTA file with annotations.
    :return: df with gene id and annotations
    """

    id= []
    annot= []
    with open(fasta, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
            annot.append(record.description)
            df = pd.DataFrame(list(zip(id, annot)), columns=["id", "desc"])
            #df["annot"] = df.apply( lambda x: x["desc"].replace(x["id"], "").strip(), axis=1) takes forever to run

    return df
dic = {}
for key, value in sp_dic.items():
    dic[key] = value

"Write annotation tables"
get_annot(dic["wb"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/wb_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["muris"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/muris_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["HIN"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/HIN_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["trepo"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/trepo_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["carpe"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/carpe_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["kbiala"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/kbiala_annot.csv",sep="\t", header="infer", index=False)
get_annot(dic["spiro"]).to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/spiro_annot.csv",sep="\t", header="infer", index=False)

"""
# I don't know why this doesn't work!
dic = {}
for key, value in sp_dic.items():
    dic[key] = get_annot(value)
    """