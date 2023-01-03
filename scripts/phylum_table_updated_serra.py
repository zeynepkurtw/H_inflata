#-Written by Begüm Serra Büyüktarakçı Nov,21
from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import argparse
import sys

def create_phylum_table(merged_fasta):

    merged  = SeqIO.parse(merged_fasta ,"fasta")
    ids = []
    for record in merged:
        if 'ANDGO' not in record.id:
            ids.append(record.id)
    print(ids)

    Entrez.email = "begumserra.buyuktarakci@icm.uu.se"  # Always tell NCBI who you are
    taxa = []
    species_name = []
    for i in ids:
        #print(i)
        handle = Entrez.efetch(db="protein", id=i, retmode="xml")
        record = Entrez.read(handle)
        taxa.append(record[0]["GBSeq_taxonomy"])
        species_name.append(record[0]["GBSeq_organism"])
    print(taxa)
    print(species_name)

    merged  = SeqIO.parse(merged_fasta ,"fasta")
    for record in merged:
        if 'ANDGO' in record.id:
            ids.append(record.id)
            taxa.append('Eukaryota; Discoba; Jakobida; Andalucina; Andaluciidae; Andalucia')
            species_name.append('Andalucia godoyi')

    df = pd.DataFrame(list(zip(ids, species_name, taxa)),
                      columns =["acc", "species", "taxa"])
    df=df.replace('Phytophthora infestans T30-4','Phytophthora infestans')
    df=df.replace('Acanthamoeba castellanii','Acanthamoeba castellanii str. Neff')
    df=df.replace('Dictyostelium discoideum','Dictyostelium discoideum AX4')

    df = pd.concat([df, df['taxa'].str.split(';', expand=True)], axis=1)
    df = df.fillna('missing')
    df = df.drop(['taxa'], axis=1)

    order = df[3]
    taxa_class = df[2]
    phylum = df[1]
    superkingdom = df[0]

    phylum_df = pd.DataFrame(list(zip(ids, species_name, order, taxa_class, phylum, superkingdom)),
                   columns =["acc", "species", "order","class","phylum", "superkingdom"])
    phylum_df.to_csv(phylum_table, sep='\t', encoding='utf-8', index=False)
    return phylum_table

###flags###
parser = argparse.ArgumentParser()
parser.add_argument("in_merged_fasta", type=str, help="merged_fasta_path")
parser.add_argument("out_phylum_table", type=str, help="phylum_table")
args = parser.parse_args()

#input
merged_fasta = args.in_merged_fasta
#output
phylum_table = args.out_phylum_table

#run
phylum_table=create_phylum_table(merged_fasta)
