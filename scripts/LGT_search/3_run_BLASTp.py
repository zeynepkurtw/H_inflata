from snakemake.shell import shell

# input,output
query = snakemake.input.query
out = snakemake.output[0]

# parameters
db_prefix = snakemake.params.db_prefix
format = snakemake.params.get("format", "")
num_threads = snakemake.params.num_threads
#max_target_seqs = snakemake.params.max_target_seqs
#max_hsps = snakemake.params.max_hsps

# command line
shell(f"""
    blastp -query {query} -db {db_prefix} -out {out} \
    -outfmt '{format}' \
    -num_threads {num_threads} \
""")
