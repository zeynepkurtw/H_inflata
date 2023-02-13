import pandas as pd
import glob
from Bio import SeqIO

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_file = "data/LGT_search/*.csv"
    fasta_file = "data/LGT_search/fasta_in/hin_trepo_cat.fasta"

    out_file = "data/LGT_search/fasta_out/{}.fasta"


# get og name from the filenames
og_dic = {}
list_files = glob.glob(og_file)
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    og_dic[i] = element


# get list of genes from OG
# get fasta from list of genes
for key, og_file in og_dic.items():

    df = pd.read_csv(og_file, sep='\t', header="infer")["id"]
    print(og_file)
    df = df.drop_duplicates().to_list()
    print("number of genes =", len(df))

    with open(fasta_file, "r") as fasta_in:
        with open(out_file.format(key), "w") as fasta_out:
            for record in SeqIO.parse(fasta_in, "fasta"):
                if record.id in df:
                    SeqIO.write(record, fasta_out, "fasta")


#concatanate fasta files










