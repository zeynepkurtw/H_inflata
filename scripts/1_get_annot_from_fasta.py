import pandas as pd
from Bio import SeqIO
import glob

path = '/Users/zeyku390/PycharmProjects/data/proteome/*.fa'
list_files = glob.glob(path)

# get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element


def get_annot(fasta):
    """
    Get a list of genes with annotations from fasta.
    :param : FASTA file with annotations.
    :return: df with gene id and annotations
    """

    id = []
    annot = []
    with open(fasta, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
            annot.append(record.description)

    df = pd.DataFrame(list(zip(id, annot)), columns=["id", "ann_f"])
    df["ann_f"] = df.apply(lambda x: x["ann_f"].replace(x["id"], "").strip(), axis=1)

    return df

dic_ann={}
for key, fasta_file in sp_dic.items():
    if key == "carpe":
        df = get_annot(fasta_file)
        df["ann_f"] = df.apply(lambda x: x["ann_f"].replace("[Carpediemonas membranifera]", "").strip(), axis=1)
        dic_ann[key] = df

    elif key == "kbiala":
        df = get_annot(fasta_file)
        df["ann_f"] = df.apply(lambda x: x["ann_f"].replace("[Kipferlia bialata]", "").strip(), axis=1)
        dic_ann[key] = df

    else:
        df = get_annot(fasta_file)

    df.to_csv(f"data/fasta_ann/{key}.csv",sep="\t", header=True, index=False)
    dic_ann[key] = df

pd.concat(dic_ann, axis=0).to_csv(f"data/fasta_ann/ann_f_cat.csv",sep="\t", header=True, index=False)


