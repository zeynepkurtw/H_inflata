from snakemake.shell import shell

input = snakemake.input
out = snakemake.output

#seq_identity = snakemake.params.seq_identity
threads = snakemake.params.threads
wildcards = snakemake.wildcards


shell(f"""cd-hit -i {input} -o {out} -c {wildcards.n} -T {threads}""")