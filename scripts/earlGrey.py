from snakemake.shell import shell

genome = snakemake.input.genome
out = snakemake.output

species = snakemake.params.species
threads = snakemake.params.threads

shell("conda activate spironucleus")
shell(f"""/Users/zeyku390/PycharmProjects/H_inflata/EarlGrey/earlGrey -g {genome} -o {out} -s {species} -t {threads}""")