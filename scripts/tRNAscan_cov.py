from snakemake.shell import shell

genome = snakemake.input.genome
tRNA = snakemake.output.tRNA
stats = snakemake.output.stats

shell(f"""tRNAscan-SE -C {genome} -o {tRNA} -m {stats}""")