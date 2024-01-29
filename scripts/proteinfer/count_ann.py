import pandas as pd
#read csv
df=pd.read_csv("resource/13_proteinfer/HIN_pred.tsv",
               sep="\t", header='infer')

#count unique values on teh first column
count = df['sequence_name'].nunique()

#print the count as number of proteins annotated by proteinfer
print(count)