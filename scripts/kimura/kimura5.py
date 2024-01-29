import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np

# Open the file
with open('resource/12_kimura_landscape/HIN.fasta.out', 'r') as f:
    lines = f.readlines()

# Remove the '*' from each line
lines = [line.rstrip(' *\n') for line in lines]

# Join the lines back into a single string
data = '\n'.join(lines)

# Load the RepeatMasker output
df = pd.read_csv(StringIO(data), sep='\s+', skiprows=3, header=None,
                 names=["SW_score", "perc_div", "perc_del", "perc_ins", "query_seq", "query_begin", "query_end", "query_left", "strand", "repeat_name", "repeat_class/family", "repeat_begin", "repeat_end", "repeat_left", "ID"])

# Calculate repeat length
df['length'] = df['query_end'] - df['query_begin'] + 1

# Convert perc_div to numeric
df['perc_div'] = pd.to_numeric(df['perc_div'])

# Remove lines with 'C' in perc_div column
df = df[df['perc_div'] != 'C']

# Remove unknown repeat elements
df = df[~df['repeat_class/family'].str.contains('Unknown')]

# Calculate Kimura divergence (CpG adjusted)
df['perc_div_cpg_adj'] = -np.log((1-df['perc_div']/100-2*df['perc_del']/100)*(1-df['perc_div']/100-2*df['perc_ins']/100)/2)*100

# Prepare data for plotting
df_grouped = df.groupby(['perc_div_cpg_adj', 'repeat_class/family']).agg({'length':'sum'}).reset_index()

# Calculate total length of the genome
genome_length = 142.7e6  # 142.7 million base pairs

# Calculate proportion
df_grouped['proportion'] = df_grouped['length'] / genome_length

# Pivot the dataframe
df_pivot = df_grouped.pivot(index='perc_div_cpg_adj', columns='repeat_class/family', values='proportion')

# Plot
plt.figure(figsize=[10,5])
df_pivot.plot(kind='bar', stacked=True)
plt.xlabel('Kimura Divergence (CpG adjusted)', fontsize=14)
plt.ylabel('Proportion of Genome', fontsize=14)
plt.title('Interspersed Repeat Landscape', fontsize=16)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()