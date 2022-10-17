from snakemake.shell import shell

input = snakemake.input
out = snakemake.output

params = snakemake.param.seq_identity
threads = snakemake.param.threads


shell(f"""cd-hit -i {input} -o {out} -c {params.seq_identity} -T {threads}""")