configfile: "env/config.yaml"

rule all:
    input:
         #"output/1_orthofinder/",
         #expand("output/2_cdhit/HIN_aa_{n}.cdhit", n=config["seq_identity"]),
         #expand("output/3_interproscan/new_sp/{n}.tsv", n=["trepo", "carpe", "kbiala"]),
         #expand("output/4_deepsig/{n}.csv", n=["HIN"])
         #expand("output/2_cdhit/{sp}_{n}.cdhit", n=config["seq_identity"], sp=config["species"])
         "/opt/zeynep/H_inflata/output/3_BLASTp/hin_trepo_cat.blastp"


rule orthofinder:
    input:
         fasta=directory("resource/1_orthofinder/fasta/"),
    output:
          directory('output/1_orthofinder/')
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/orthofinder.py"

rule orthofinder_rerun:
    input:
         new_sp=directory("resource/1_orthofinder/new_sp/"),
         #old_sp= directory("resource/1_orthofinder/WorkingDirectory_5sp/OrthoFinder/Results_Oct17_2/WorkingDirectory_7p/")
         old_sp=directory("resource/1_orthofinder/WorkingDirectory/OrthoFinder/Results_Oct17_2/WorkingDirectory/")
    output:
          directory('output/1_orthofinder/new_sp')
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/orthofinder_rerun.py"

rule cdhit:
    input: "resource/2_cdhit/{sp}_aa.fasta"
    params:
          threads=8
    output: "output/2_cdhit/{sp}_{n}.cdhit"
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/cdhit.py"

rule interproscan:
    input:
         "resource/3_interproscan/new_sp/{sp}_aa.fasta"
    params:
          threads=32,
          output_folder="output/3_interproscan/new_sp"
    output:
          "output/3_interproscan/new_sp/{sp}.tsv"
    #directory("output/3_interproscan/new_sp/{sp}"),
    script: "scripts/interproscan.py"

rule trepo_list:
    input:
         og="output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt"
    output:
          "/Users/zeyku390/PycharmProjects/H_inflata/plots/upset_trepo.png"
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/10_LGT_upset.py"

rule blastp:
    input:
            query = "/opt/zeynep/H_inflata/resource/6_BLASTp/hin_trepo_cat.fasta",
            db = "/data/zeynep/databases"
    output:
            "/opt/zeynep/H_inflata/output/3_BLASTp/hin_trepo_cat.blastp"
    params:
          format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids",
          num_threads=30,
          evalue=1e-10,
          perc_identity= "perc_identity",
          db_prefix="/data/zeynep/databases/nr"
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/LGT_search/3_run_BLASTp.py"

"""
rule deepsig:
    input: "resource/4_deepsig/{sp}_aa.fasta"
    params: threads= 30
    output: "output/4_deepsig/{sp}.csv"
    conda: "env/hinflata.yaml"
    script: "scripts/deepsig.py"


rule template:
    input: "resource/"
    params: ""
    output: "output/"
    conda: "env/hinflata.yaml"
    script: ""
"""
