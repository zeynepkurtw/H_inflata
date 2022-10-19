from snakemake.shell import shell

input = snakemake.input
output = snakemake.output

#seq_identity = snakemake.params.seq_identity
threads = snakemake.params.threads


shell(f"""deepsig.py -f {input} -o {output} -k euk -a {threads}""")