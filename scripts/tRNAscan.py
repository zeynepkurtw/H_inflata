from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
stats = snakemake.output.stats

#seq_identity = snakemake.params.seq_identity
threads = snakemake.params.threads


shell(f"""tRNAscan-se {input} -o {output} -m {stats}""")