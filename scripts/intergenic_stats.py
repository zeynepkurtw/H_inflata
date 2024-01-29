import pandas as pd

import pandas as pd

def calculate_lengths(intergenic_bed, genes):
    intergenic_df = pd.read_csv(intergenic_bed, sep='\t', header=None)
    genes_df = pd.read_csv(genes, sep='\t', header=None)

    intergenic_lengths = intergenic_df[2] - intergenic_df[1] + 1
    gene_lengths = genes_df[4] - genes_df[3] + 1

    mean_intergenic_length = intergenic_lengths.mean()
    median_intergenic_length = intergenic_lengths.median()
    mean_protein_length = (gene_lengths / 3).mean()
    median_protein_length = (gene_lengths / 3).median()
    total_intergenic_length = intergenic_lengths.sum()
    gene_density = gene_lengths.sum() / (gene_lengths.sum() + intergenic_lengths.sum())
    num_protein_coding_genes = len(genes_df)

    return (mean_intergenic_length, median_intergenic_length,
            mean_protein_length, median_protein_length,
            total_intergenic_length, gene_density, num_protein_coding_genes,
            intergenic_df, genes_df)  # Return intergenic_df and genes_df

sp_list = ["HIN", "carpe", "muris", "spiro", "kbiala"]

for sp in sp_list:
    intergenic_bed = f'output/6_bedtools/{sp}.intergenic.bed'
    genes = f'output/6_bedtools/{sp}.CDS.sorted.gff'

    results = calculate_lengths(intergenic_bed, genes)

    mean_intergenic, median_intergenic, mean_protein, median_protein, total_intergenic, gene_density, num_protein_coding_genes, intergenic_df, genes_df = results

    print(f'For {sp}:')
    print(f'Number of protein coding genes: {num_protein_coding_genes}')
    print(f'Mean intergenic length: {mean_intergenic}')
    print(f'Median intergenic length: {median_intergenic}')
    print(f'Mean protein length: {mean_protein}')
    print(f'Median protein length: {median_protein}')
    print(f'Total intergenic length: {total_intergenic}')
    print(f'Gene density: {gene_density}\n')



    """
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
    
    # Calculate gene density as the ratio of the total length of genes
    # to the total length of the genome
    genome_length = gene_lengths.sum() + intergenic_lengths.sum()
    gene_density = gene_lengths.sum() / genome_length
    print(f'Gene density: {gene_density}')


    """


