from snakemake.shell import shell

index_genome = snakemake.input.index_genome
genes = snakemake.input.genes

cut_genome = snakemake.output.cut_genome
genes_sorted = snakemake.output.genes_sorted
intergenic= snakemake.output.intergenic

# First, two columns from index genome
# Second, sort genes.bed, in the subshell, first field contig order, second field start positions numerically
# Third, complement to extract intergenic regions with given genes.gff (genes.bed)

shell(f"cut -f1,2,3 {index_genome} | sort -V -k1,1 > {cut_genome}")
shell(f"sort -V -k1,1 -k4,4n -k5,5n {genes} > {genes_sorted}")
#shell(f"bedtools sort -i {genes} > {genes_sorted}")
shell(f"bedtools complement -i {genes_sorted} -g {cut_genome} > {intergenic}")
