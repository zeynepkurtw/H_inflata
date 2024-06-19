import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Setup
sns.set_style("whitegrid", {'axes.grid' : False})
cont_blast = "resource/15_assembly_cont/prok_cont_hybrid_polished_masurca.blastn"

# Load and process the data
data = pd.read_csv(cont_blast, header=None, sep='\t').sort_values(by=[3], ascending=False)
df = data.iloc[:, [0, 2, 3, 6, 7, 10, 12, 14]]
df = df.rename(columns={0: "qseqid", 2: "pident", 3: "length", 6: "qstart", 7: "qend", 10: "e_value", 12: "qlen", 14: "stitle"})
df['stitle_short'] = df['stitle'].str.split(" ", expand=True)[[0, 1]].agg(' '.join, axis=1)
df = df.drop(columns="stitle")

# Filter the data
df_filtered = df[df.e_value < pow(10, -10)]
df_cleaned = df_filtered[~df_filtered.qseqid.isin(["contig_573", "contig_427", "contig_758"])]

# Create a figure with two subplots (side by side)
fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
sns.scatterplot(ax=axs[0], data=df_filtered, x="qlen", y="length", hue="pident").set_title("Before Cleaning")
sns.scatterplot(ax=axs[1], data=df_cleaned, x="qlen", y="length", hue="pident", s=100).set_title("After Cleaning")

# Adjust legend and axis formatting
axs[0].get_legend().remove()  # Remove the legend from the first plot
axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, title='pi')
# Rotate x-axis labels for both subplots to 45 degrees
plt.setp(axs[0].get_xticklabels(), rotation=45)
plt.setp(axs[1].get_xticklabels(), rotation=45)
  # Rotate x-axis labels for both subplots

# Adjust labels
axs[0].set_xlabel("Contig")
axs[0].set_ylabel("Alignment")
axs[1].set_xlabel("Contig")
axs[1].set_ylabel("Alignment")

# Adjust layout and save
plt.tight_layout()
plt.savefig("plots/assembly_cont/combined_plots.png", dpi=300, bbox_inches='tight')
plt.show()


