from snakemake.shell import shell

input = snakemake.input
out = snakemake.params.output_folder

threads = snakemake.params.threads


shell(f"""/data/zeynep/interproscan-5.47-82.0/interproscan.sh -i {input} -d {out} -f gff3,tsv,json -iprlookup -goterms --pathways -cpu {threads} """)