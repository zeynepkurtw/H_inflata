import pandas as pd

barrnap ="/Users/zeyku390/PycharmProjects/H_inflata/output/14_barRNAp/HIN.rrna.gff"
#barrnap ="/Users/zeyku390/PycharmProjects/H_inflata/output/14_barRNAp/muris.rrna.gff"


df = pd.read_csv(barrnap, sep='\t', comment='#', header=None)
df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
df['length'] = df['end'] - df['start']
total_length = df['length'].sum()
num_rrna = df.shape[0]
average_length = total_length / num_rrna

#percantage of total length of rRNA in the genome
total_genome_length = 142.6e6 #HIN
#total_genome_length = 9.8e6 #muris
percentage = total_length / total_genome_length * 100
print (f"Percentage of total length of rRNA in the genome: {percentage:.2f}%")




print(f"Total length of rRNA: {total_length}")
print(f"Number of rRNA: {num_rrna}")
