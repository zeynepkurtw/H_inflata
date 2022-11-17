configfile: "env/config.yaml"

rule all:
    input:
        #"output/1_orthofinder",
        expand("output/2_cdhit/HIN_aa_{n}.cdhit", n=config["seq_identity"]),
        #expand("output/3_interproscan/new_sp/{n}.tsv", n=["trepo", "carpe", "kbiala"]),
        #expand("output/4_deepsig/{n}.csv", n=["HIN"])

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
        threads= 8
    output: "output/2_cdhit/{hin}_{n}.cdhit"
    conda:
        "env/hinflata.yaml"
    script:
        "scripts/cdhit.py"

rule interproscan:
    input:
        "resource/3_interproscan/new_sp/{sp}_aa.fasta"
    params:
        threads= 32 ,
        output_folder = "output/3_interproscan/new_sp"
    output:
        "output/3_interproscan/new_sp/{sp}.tsv"
        #directory("output/3_interproscan/new_sp/{sp}"),

    script: "scripts/interproscan.py"

rule deepsig:
    input: "resource/4_deepsig/{sp}_aa.fasta"
    params: threads= 30
    output: "output/4_deepsig/{sp}.csv"
    conda: "env/hinflata.yaml"
    script: "scripts/deepsig.py"

"""
rule template:
    input: "resource/"
    params: ""
    output: "output/"
    conda: "env/hinflata.yaml"
    script: ""
"""