from snakemake.shell import shell

input = snakemake.input
output_folder = snakemake.params.output_folder

threads = snakemake.params.threads

#sed 's/*//' trepo.faa > trepo_aa.fasta
shell(f"""/data/zeynep/interproscan-5.47-82.0/intercd ../.proscan.sh -i {input} -d {output_folder} -f gff3,tsv,json -iprlookup -goterms --pathways -cpu {threads} """)