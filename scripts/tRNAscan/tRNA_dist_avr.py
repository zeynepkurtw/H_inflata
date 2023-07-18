import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("output/5_tRNAscan/HIN.tRNAscan", sep="\t", header=None, skiprows=3, usecols=[0,1,2,3,4,5])

#df.columns = df.columns.str.strip()
df.columns =["Name","tRNA", "Begin", "End","Type","Codon"]
df= df.drop_duplicates()

# Count the number of tRNAs for each contig
tRNA_counts = df['Name'].value_counts()

# Compute the average number of tRNAs
average_tRNA = tRNA_counts.mean()

# Find the contig that is closest to the average
contig_with_average_tRNAs = (tRNA_counts - average_tRNA).abs().idxmin()




# Filter the dataframe to only include rows from the selected contig
df_filtered = df[df['Name'] == contig_with_average_tRNAs]
df_filtered['Begin'] = pd.to_numeric(df_filtered['Begin'], errors='coerce')
df_filtered['End'] = pd.to_numeric(df_filtered['End'], errors='coerce')
df_filtered = df_filtered.sort_values(by="Begin", ascending=True)

# Generate y values for each tRNA
df_filtered['y'] = range(1, len(df_filtered) + 1)

# Create a new figure
plt.figure(figsize=(10, 6))

# Plot the Begin and End positions
plt.scatter(df_filtered['Begin'], df_filtered['y'], alpha=0.5, color='blue', label='Begin')
plt.scatter(df_filtered['End'], df_filtered['y'], alpha=0.5, color='red', label='End')

# Draw lines between Begin and End positions
for _, row in df_filtered.iterrows():
    plt.plot([row['Begin'], row['End']], [row['y'], row['y']], color='black', linewidth=0.5)

plt.xlabel('Position')
plt.ylabel('tRNA')
plt.legend()
plt.title('tRNA Begin and End Positions')

plt.savefig(out, format="svg", bbox_inches='tight', dpi=300)

plt.show()
