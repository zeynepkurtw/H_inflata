configfile: "env/config.yaml"

rule all:
    input:
        #"output/1_orthofinder",
        expand("output/2_cdhit/HIN_{n}.cdhit", n=config["seq_identity"])

rule orthofinder:
    input:
        new_sp= directory("resource/1_orthofinder/new_sp"),
        old_sp= directory("resource/1_orthofinder/WorkingDirectory")
    output:
        directory('output/1_orthofinder')
    conda:
        "env/hinflata.yaml"
    script:
        "scripts/orthofinder.py"

rule cdhit:
    input:"resource/2_cdhit/{hin}.fasta"
    params:
        threads= 30
    output: "output/2_cdhit/{hin}_{n}.cdhit"
    conda:
        "env/hinflata.yaml"
    script:
        "scripts/cdhit.py"