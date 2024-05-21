import pandas as pd
#read csv
interpro ="resource/3_interproscan/HIN.tsv",


df = pd.read_csv(interpro, sep="\t", header=None, names=list(range(0, 15)),
                           engine='python', quoting=3)[[0, 3, 4, 5, 11, 12]]
df= df.dropna().drop_duplicates().rename(
    columns={0: "id", 3: "db", 4: "db_acc", 5: "ann_db", 11: "ipr", 12: "ann_inter"})

#count unique values on teh first column
id = len(df['id'].drop_duplicates())
#count PFAM domains
pfam = df['db_acc'].str.contains("PF").sum()
#count GO terms

LRR = df['description'].str.contains("LRR|leucine").sum()
cysteine = df['description'].str.contains("cysteine").sum()

#count PFAM with confiedence level more tahn 95
pfam95 = df[(df['predicted_label'].str.contains("Pfam")) & (df['confidence'] > 0.95)].shape[0]

#update the script to count InterPro annotations



#print the count as number of proteins annotated by proteinfer
print("id", id)
print("PFAM", pfam)
print("GO", go)
print("LRR", LRR)
print("PFAM95", pfam95)
print("cysteine", cysteine)