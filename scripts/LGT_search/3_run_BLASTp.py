from snakemake.shell import shell

# input,output
query = snakemake.input.query
out = snakemake.output[0]

# parameters
db_prefix = snakemake.params.db_prefix
perc_identity = snakemake.params.perc_identity
format = snakemake.params.format
num_threads = snakemake.params.num_threads

# command line
shell(f"""
    blastp -query {query} -db {db_prefix} -out {out} \
    -perc_identity {perc_identity} \
    -outfmt {format} \
    -num_threads {num_threads}
""")