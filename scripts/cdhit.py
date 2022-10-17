from snakemake.shell import shell

input = snakemake.input
out = snakemake.output

seq_identity = snakemake.params.seq_identity
threads = snakemake.params.threads


shell(f"""cd-hit -i {input} -o {out} -c {seq_identity} -T {threads}""")