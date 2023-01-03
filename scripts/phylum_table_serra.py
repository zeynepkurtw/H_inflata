#-Written by Begüm Serra Büyüktarakçı Nov,21
from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import argparse
import sys

def extract_ids(blast_file):
    """
    extract protein ids from blast file (for outgroups)
    """
    blast=pd.read_csv(blast_file, sep='\t', header=None)
    blast_col = blast[1].drop_duplicates()
    ids_outgroups = blast_col
    return ids_outgroups

def write_ids(ids_outgroups):
    """
    write protein ids (outgroups) into o file
    """
    ids_outgroups_list=[]
    for i in ids_outgroups:
        ids_outgroups_list.append(i)
    #print(ids_outgroups_list)
    textfile = open(ids_outgroup_file, "w")
    for element in ids_outgroups_list:
        textfile.write(element + "\n")
    textfile.close()
    return ids_outgroups_list

def extract_taxa(ids_outgroups):
    """
    extract taxa from ncbi for both outgroups and orthogroups
    """
    ids_outgroups_list = ids_outgroups.values.tolist()
    fasta_ortho = SeqIO.parse(fasta_ortho_file,"fasta")
    ids_orthogroups = []
    for record in fasta_ortho:
        ids_orthogroups.append(record.id)
        #ids = ids_orthogroups + ids_outgroups_list
        ids = ids_orthogroups

    #print(ids)
    ids_df = pd.DataFrame(ids)
    ids_df.to_csv(ids_outgroup_orthogroup, sep='\t', encoding='utf-8', index=False)

    Entrez.email = "begumserra.buyuktarakci@icm.uu.se"  # Always tell NCBI who you are
    taxa = []
    species_name = []
    for i in ids:
        #print(i)
        handle = Entrez.efetch(db="protein", id=i, retmode="xml")
        record = Entrez.read(handle)
        taxa.append(record[0]["GBSeq_taxonomy"])
        species_name.append(record[0]["GBSeq_organism"])
        #print(species_name)
        #print(type(species_name))
        taxa_df = pd.DataFrame(taxa)
    taxa_df.to_csv(taxa_file, sep='\t', encoding='utf-8', index=False, header=None)
    return (taxa_file, ids_outgroup_orthogroup, species_name)

def create_phylum_table(taxa_file, ids_outgroup_orthogroup, species_name):
    """
    create phylum table with protein ids, species, phylum and superkingdom
    """
    protein_ids = pd.read_csv(ids_outgroup_orthogroup, sep='\t', encoding='utf-8')
    protein_ids = protein_ids["0"]
    column_title=["col0","col1","col2","col3","col4","col5","col6","col7","col8","col9","col10","col11","col12","col13","col14","col115","col16"]
    taxa=pd.read_csv(taxa_file, sep=";", header=None, names=column_title)
    order = taxa["col3"]
    taxa_class = taxa["col2"]
    phylum = taxa["col1"]
    superkingdom = taxa["col0"]
    print(phylum)
    print(superkingdom)

    phylum_df = pd.DataFrame(list(zip(protein_ids, species_name, order, taxa_class, phylum, superkingdom)),
                   columns =["acc", "species", "order","class","phylum", "superkingdom"])
    phylum_df.to_csv(phylum_table, sep='\t', encoding='utf-8', index=False)
    return phylum_table

###flags###
parser = argparse.ArgumentParser()
parser.add_argument("in_blast_file", type=str, help="blast_file_path")
parser.add_argument("in_fasta_ortho_file", type=str, help="fasta_ortho_file_path")
parser.add_argument("in_seqs_outgroup_file", type=str, help="seqs_outgroup_file_path")

parser.add_argument("out_ids_outgroup_file", type=str, help="ids_outgroup_file_path")
parser.add_argument("out_taxa_file", type=str, help="taxa_file")
parser.add_argument("out_phylum_table", type=str, help="phylum_table")
parser.add_argument("out_ids_outgroup_orthogroup", type=str, help="ids_outgroup_orthogroup")

args = parser.parse_args()
#input
blast_file = args.in_blast_file
fasta_ortho_file = args.in_fasta_ortho_file
seqs_outgroup_file = args.in_seqs_outgroup_file

#output
ids_outgroup_file = args.out_ids_outgroup_file
taxa_file= args.out_taxa_file
phylum_table = args.out_phylum_table
ids_outgroup_orthogroup = args.out_ids_outgroup_orthogroup


#run
ids_outgroups = extract_ids(blast_file)
ids_outgroup_file = write_ids(ids_outgroups)
taxa_file, ids_outgroup_orthogroup, species_name=extract_taxa(ids_outgroups)
phylum_table=create_phylum_table(taxa_file, ids_outgroup_orthogroup, species_name)
