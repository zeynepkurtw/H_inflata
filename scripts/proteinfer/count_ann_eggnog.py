import pandas as pd
#read csv
df=pd.read_csv("resource/13_proteinfer/spiro_pred.tsv",
               sep="\t", header='infer')

#count unique values on teh first column
id = len(df['sequence_name'].drop_duplicates())
#count PFAM domains
pfam = df['predicted_label'].str.contains("Pfam").sum()
#count GO terms
go = df['predicted_label'].str.contains("GO").sum()
LRR = df['description'].str.contains("LRR|leucine").sum()
cysteine = df['description'].str.contains("cysteine").sum()

#count PFAM with confiedence level more tahn 95
pfam95 = df[(df['predicted_label'].str.contains("Pfam")) & (df['confidence'] > 0.95)].shape[0]

#print the count as number of proteins annotated by proteinfer
print("id", id)
print("PFAM", pfam)
print("GO", go)
print("LRR", LRR)
print("PFAM95", pfam95)
print("cysteine", cysteine)