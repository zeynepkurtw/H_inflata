from snakemake.shell import shell

index_genome = snakemake.input.index_genome
genes = snakemake.input.genes

intergenic= snakemake.output.intergenic


shell(f"bedtools complement -i {genes} -g {index_genome} > {intergenic}")
