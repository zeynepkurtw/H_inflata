from snakemake.shell import shell

genome = snakemake.input.genome
tRNA = snakemake.output.tRNA
stats = snakemake.output.stats

shell(f"""trnascan-1.4 {genome} -o {tRNA} -m {stats}""")