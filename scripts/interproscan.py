from snakemake.shell import shell

input = snakemake.input
output = snakemake.output

threads = snakemake.params.threads


shell(f"""/data/zeynep/interproscan-5.47-82.0/interproscan.sh -i {input} -d {output} -f gff3,tsv,json -iprlookup -goterms --pathways -cpu {threads} """)