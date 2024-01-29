import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Or try 'TkAgg', 'Qt5Agg', etc.
import matplotlib.pyplot as plt

input= "output/5_tRNAscan/HIN.tRNAscan"
out1= "plots/tRNAscan/HIN_tRNA_all.svg"
out2= "plots/tRNAscan/HIN_tRNA_all.png"
out3= "plots/tRNAscan/HIN_tRNA_all.tiff"

# Assuming you have a DataFrame df, with 'start' and 'end' columns for tRNA locations
df = pd.read_csv(input, sep="\t", header=None, skiprows=3, usecols=[0,1,2,3,4,5])

df.columns =["Name","tRNA", "Begin", "End","Type","Codon"]
df= df.drop_duplicates()

df['Begin'] = pd.to_numeric(df['Begin'], errors='coerce')
df['End'] = pd.to_numeric(df['End'], errors='coerce')

# Generate an arbitrary y value for each tRNA for visualization purposes
df['y'] = range(1, len(df) + 1)

# Create the scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(df['Begin'], df['y'], alpha=0.5, color="#7A93A2", label='Begin')
#plt.scatter(df['End'], df['y'], alpha=0.5, color='red', label='End')


plt.xlabel('Position in Genome')
plt.ylabel('tRNA')
plt.title('Distribution of tRNAs Across the Genome')

plt.savefig(out1, format="svg", bbox_inches='tight', dpi=300)
plt.savefig(out2, format="png", bbox_inches='tight', dpi=300)
plt.savefig(out3, format="tiff", bbox_inches='tight', dpi=300)
plt.show()

