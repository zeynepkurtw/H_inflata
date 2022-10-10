from snakemake.shell import shell

# input,output
new_sp = snakemake.input.new_sp
old_sp =snakemake.input
out = snakemake.output



#OrthoFinder allows you to add extra species without re-running the previously computed BLAST searches:
shell(f"""orthofinder -b {old_sp} -f {new_sp} -S blast""")