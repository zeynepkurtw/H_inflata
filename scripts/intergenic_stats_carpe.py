import pandas as pd

intergenic_bed = 'output/6_bedtools/carpe.intergenic.bed'
genes = "output/6_bedtools/carpe.CDS.sorted.gff"

# Load the intergenic.bed and genes.bed file
intergenic_df = pd.read_csv(intergenic_bed, sep='\t', header=None)
genes_df = pd.read_csv(genes, sep='\t', header=None)

# Columns 1 and 2 contain the start and end of the regions
# So, we subtract the start from the end to get the lengths
intergenic_lengths = intergenic_df[2] - intergenic_df[1] + 1
gene_lengths = genes_df[4] - genes_df[3] + 1

# Calculate the mean and median of the intergenic region lengths
mean_intergenic_length = intergenic_lengths.mean()
median_intergenic_length = intergenic_lengths.median()
print(f'Mean intergenic length: {mean_intergenic_length}')
print(f'Median intergenic length: {median_intergenic_length}')

# Calculate the mean and median of the protein lengths
# Assuming that 3 nucleotides code for one amino acid
protein_lengths = gene_lengths / 3
mean_protein_length = protein_lengths.mean()
median_protein_length = protein_lengths.median()
print(f'Mean protein length: {mean_protein_length}')
print(f'Median protein length: {median_protein_length}')

# Calculate the total intergenic length
total_intergenic_length = intergenic_lengths.sum()
print(f'Total intergenic length: {total_intergenic_length}')

# Calculate the average intergenic length
average_intergenic_length = intergenic_lengths.mean()
print(f'Average intergenic length: {average_intergenic_length}')

# Calculate gene density as the ratio of the total length of genes
# to the total length of the genome
genome_length = gene_lengths.sum() + intergenic_lengths.sum()
gene_density = gene_lengths.sum() / genome_length
print(f'Gene density: {gene_density}')