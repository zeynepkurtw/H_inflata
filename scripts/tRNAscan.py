from snakemake.shell import shell

genome = snakemake.input.genome
tRNA = snakemake.output.tRNA
stats = snakemake.output.stats
gff= snakemake.output.gff

#shell(f"""tRNAscan-SE {genome} -o {tRNA} -m {stats}""")
shell(f"""tRNAscan-SE {genome} -o {tRNA} -f {gff} -m {stats}""")

