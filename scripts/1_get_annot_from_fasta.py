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
    :param old: FASTA file with annotations.
    :return: df with gene id and annotations
    """

    id = []
    annot = []
    with open(fasta, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            id.append(record.id)
            annot.append(record.description)

    df = pd.DataFrame(list(zip(id, annot)), columns=["id", "desc"])
    df["desc"] = df.apply(lambda x: x["desc"].replace(x["id"], "").strip(), axis=1)

    return df


for key, fasta_file in sp_dic.items():
    get_annot(fasta_file).to_csv(f"/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/{key}_annot.csv",
                                 sep="\t", header=True, index=False)

