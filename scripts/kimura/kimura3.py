import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

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

# Prepare data for plotting
df_grouped = df.groupby(['perc_div']).agg({'length':'sum'}).reset_index()

# Plot
plt.figure(figsize=[10,5])
plt.bar(df_grouped['perc_div'], df_grouped['length'], color='blue', edgecolor='green')
plt.xlabel('Kimura Divergence', fontsize=14)
plt.ylabel('Total Repeat Length (bp)', fontsize=14)
plt.title('Kimura Divergence Landscape', fontsize=16)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
