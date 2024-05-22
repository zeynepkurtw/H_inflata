import pandas as pd

superfamily= "LRR"

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    interproscan_file = "jupyter/Annas_script/multifamily/input/HIN.tsv"

    interpro_pfam = "scripts/paper3/pfam_barplot/interpro_HIN_LRR_pfam.csv"

interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession',
                           'Signature description','Start location', 'Stop location', 'Score', 'Status',
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
                           'optional: GO terms', 'optional: Pathway annotation']

df = pd.read_csv(interproscan_file, sep='\t', header=None, names = interproscan_header)

#get all proteins with IPR032675= LRR superfamily

df_=df[df["Interpro annotations - accession"] == "IPR032675"]
df_.groupby("Interpro annotations - accession").count()


#df_pfam.to_csv(interpro_pfam, index=False, sep="\t")


print ( "number of ", superfamily, "protein = ",  len(df_.drop_duplicates("protein_id")))
print ( "number of ", superfamily, "domain = ",  len(df_.drop_duplicates()))




