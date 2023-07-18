import pandas as pd
# claculate mean-median protein size
# integenic regions
# introns

#read gff or protein multifasta file

input="resource/9_genome_stats/HIN_single-exon_79332.gff"

out= "data/genome_stats"


# Load the GFF file
df = pd.read_csv(input, sep='\t', comment='#', header=None,
                 names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])[:100]

# Filter the dataframe to only include gene features
genes = df[df['type'] == 'CDS'].copy()

# Calculate the size of each gene
genes['size'] = genes['end'] - genes['start'] + 1

# Calculate protein sizes in amino acids
genes['protein_size'] = genes['size'] // 3

# Calculate mean protein size
mean_protein_size = genes['protein_size'].mean()

# Calculate median protein size
median_protein_size = genes['protein_size'].median()

# Calculate total gene size
total_gene_size = genes['size'].sum()

# Create a sorted list of all gene and intergenic region boundaries
boundaries = sorted(genes['start'].tolist() + genes['end'].tolist())

# Calculate the sizes of intergenic regions
intergenic_sizes = [boundaries[i+1] - boundaries[i] - 1 for i in range(0, len(boundaries)-1, 2)]

# Calculate total size of intergenic regions
total_intergenic_region = sum(intergenic_sizes)

# Calculate mean intergenic size
mean_intergenic_size = total_intergenic_region / len(intergenic_sizes) if intergenic_sizes else 0

# Calculate median intergenic size
intergenic_sizes.sort()
middle = len(intergenic_sizes) // 2
if len(intergenic_sizes) % 2 == 0:
    median_intergenic_size = (intergenic_sizes[middle - 1] + intergenic_sizes[middle]) / 2
else:
    median_intergenic_size = intergenic_sizes[middle]

# Calculate coding density
genome_size = 143 * (10 ** 6) # 143 Mbp
coding_density = total_gene_size / genome_size

print(f"Mean protein size: {mean_protein_size} amino acids")
print(f"Median protein size: {median_protein_size} amino acids")
print(f"Total size of intergenic regions: {total_intergenic_region} base pairs")
print(f"Mean intergenic size: {mean_intergenic_size} base pairs")
print(f"Median intergenic size: {median_intergenic_size} base pairs")
print(f"Coding density: {coding_density}")