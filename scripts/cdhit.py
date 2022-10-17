from snakemake.shell import shell

input = snakemake.input
out = snakemake.output[0]

param = snakemake.param.seq_identity
threads = snakemake.param.threads


shell(f"""cd-hit -i {input} -o {out} -c {param.seq_identity} -T {threads}""")