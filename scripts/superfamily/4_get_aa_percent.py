import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import numpy as np
from collections import defaultdict

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    path = '/Users/zeyku390/PycharmProjects/data/proteome/*.fa'
    list_files = glob.glob(path)
    out_file = "data/superfamily/aa_percent-{}.csv"

# get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element


def get_seq_dict(fasta):
    ''' Return sequence dictionary '''

    seq_d = defaultdict(Seq)
    with open(fasta, 'r') as fasta_in:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            seq_d[record.id] = record.seq

    return seq_d


def get_seq_stats(seq_d):
    ''' Return list of aa count and percentage, and sequences of given sequences'''

    lens = {}
    aa_L = {}
    aa_C = {}
    p_L = {}
    p_C = {}

    for seqid, seq in seq_d.items():
        aa_L[seqid] = seq.count('L')
        aa_C[seqid] = seq.count('C')
        lens[seqid] = len(seq)
        p_L[seqid] = (seq.count('L') / len(seq)) * 100
        p_C[seqid] = (seq.count('C') / len(seq)) * 100

    df = pd.DataFrame.from_dict([p_L, p_C, lens]).T.reset_index()
    df = df.rename(columns={"index": "id", 0: "perc_L", 1: "perc_C", 2: "len"})
    df = df.sort_values(by=['perc_L', "perc_C"], ascending=False)
    #df[["perc_L", "perc_C"]] = df[["perc_L", "perc_C"]].apply(lambda x: round(x, 2))
    df[["perc_L", "perc_C"]] = df[["perc_L", "perc_C"]].round(2)
    df[["len"]] = df[["len"]].apply(pd.to_numeric, downcast='integer')

    return df


for key, value in sp_dic.items():
    seq_d = get_seq_dict(value)
    sp_dic[key] = get_seq_stats(seq_d)
    sp_dic[key].to_csv(out_file.format(key), index=False, sep="\t")

    for aa in ['Leucine', 'Cysteine']:
        print(key + f" mean {aa} = ", round(sp_dic[key][f"perc_{aa[0]}"].mean(), 2))
        print(key + f" median {aa} = ", round(sp_dic[key][f"perc_{aa[0]}"].median(), 2))
        print(key + f" max {aa} = ", round(sp_dic[key][f"perc_{aa[0]}"].max(), 2))

    print("")
