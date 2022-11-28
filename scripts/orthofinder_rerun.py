from snakemake.shell import shell

# input,output
new_sp = snakemake.input.new_sp
old_sp =snakemake.input.old_sp
out = snakemake.output



#OrthoFinder allows you to add extra species without re-running the previously computed BLAST searches:
shell(f"""orthofinder -f {new_sp} -b {old_sp} -S blast""")
