from snakemake.shell import shell

# input,output
fasta = snakemake.input.new_sp
out = snakemake.output



# -og : stop after inferring orthogroups
shell(f"""orthofinder -f {fasta} -og -S blast""")
