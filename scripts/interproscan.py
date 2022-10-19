from snakemake.shell import shell

input = snakemake.input
output_folder = snakemake.params.output_folder

threads = snakemake.params.threads

shell(f"""sed 's/*//' resource/3_interproscan/new_sp/trepo_aa.fasta > resource/3_interproscan/new_sp/trepo_aa.fasta """)
shell(f"""/data/zeynep/interproscan-5.47-82.0/interproscan.sh -i {input} -d {output_folder} -f gff3,tsv,json -iprlookup -goterms --pathways -cpu {threads} """)