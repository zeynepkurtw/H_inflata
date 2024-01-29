from snakemake.shell import shell

genome = snakemake.input.genome
index_genome = snakemake.output.index_genome

shell(f"samtools faidx {genome} -o {index_genome}")
