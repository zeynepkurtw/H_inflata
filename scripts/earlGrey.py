from snakemake.shell import shell

genome = snakemake.input.genome
out = snakemake.output.out

species = snakemake.params.species
threads = snakemake.params.threads

shell(f"""earlGrey -g {genome} -o {out} -s {species} -t {threads}""")