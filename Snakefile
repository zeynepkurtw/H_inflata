configfile: "config.yaml"

rule all:
    input:
        "output/1_orthofinder/OrthoFinder"


rule orthofinder:
    input:
        new_sp= "resource/1_orthofinder/new_sp/carpe.faa",
        old_sp= directory("resource/1_orthofinder/WorkingDirectory")
    output:
        directory('output/1_orthofinder/OrthoFinder')
    conda:
        "env/hinflata.yaml"
    script:
        "scripts/orthofinder.py"