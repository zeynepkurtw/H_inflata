import pandas as pd
import matplotlib.pyplot as plt

input= "output/5_tRNAscan/HIN.tRNAscan"
out= "plots/tRNAscan/HIN_tRNA_max.svg"

# Assuming you have a DataFrame df, with 'contig', 'start', and 'end' columns
df = pd.read_csv(input, sep="\t", header=None, skiprows=3, usecols=[0,1,2,3,4,5])

#df.columns = df.columns.str.strip()
df.columns =["Name","tRNA", "Begin", "End","Type","Codon"]
df= df.drop_duplicates()

# Group by contig and get the one with maximum count
most_abundant_contig = df['Name'].value_counts().idxmax()

# Filter the dataframe to only include rows from the most abundant contig
df_filtered = df[df['Name'] == most_abundant_contig]
df_filtered['Begin'] = pd.to_numeric(df_filtered['Begin'], errors='coerce')
df_filtered['End'] = pd.to_numeric(df_filtered['End'], errors='coerce')

#df_filtered = df_filtered.sort_values(by="Begin", ascending=True)

# Generate an arbitrary y value for each tRNA for visualization purposes
df_filtered['y'] = range(1, len(df_filtered) + 1)

# Create the scatter plot
plt.figure(figsize=(20, 10))

# Plot the Begin and End positions
plt.scatter(df_filtered['Begin'], df_filtered['y'], alpha=0.5, color='blue', label='Begin')
#plt.scatter(df_filtered['End'], df_filtered['y'], alpha=0.5, color='red', label='End')


# Draw lines between Begin and End positions
for _, row in df_filtered.iterrows():
    plt.plot([row['Begin'], row['End']], [row['y'], row['y']], color='black', linewidth=0.5)


plt.xlabel('Position in Contig')
plt.xticks(rotation=45)
plt.ylabel('tRNA')
plt.title(f'Distribution of tRNAs in Contig: {most_abundant_contig}')

plt.savefig(out, format="svg", bbox_inches='tight', dpi=300)
plt.show()
