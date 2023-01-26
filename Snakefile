configfile: "env/config.yaml"
n_partitions_blastp = config["n_partitions_blastp"]


rule all:
    input:
         #"output/1_orthofinder/",
         #expand("output/2_cdhit/HIN_aa_{n}.cdhit", n=config["seq_identity"]),
         #expand("output/3_interproscan/new_sp/{n}.tsv", n=["trepo", "carpe", "kbiala"]),
         #expand("output/4_deepsig/{n}.csv", n=["HIN"])
         #expand("output/2_cdhit/{sp}_{n}.cdhit", n=config["seq_identity"], sp=config["species"])
         #"/opt/zeynep/H_inflata/output/3_BLASTp/hin_trepo_cat.blastp"
         expand("output/3_BLASTp/{file}_{i_partition}.blastp",file=["ss_trepo"], i_partition=range(n_partitions_blastp))


rule orthofinder:
    input:
         fasta="resource/1_orthofinder/fasta/",
    output:
          directory('output/1_orthofinder/')
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/orthofinder.py"

rule orthofinder_rerun:
    input:
         new_sp="resource/1_orthofinder/new_sp/",
         #old_sp="resource/1_orthofinder/WorkingDirectory_5sp/OrthoFinder/Results_Oct17_2/WorkingDirectory_7p/"
         old_sp="resource/1_orthofinder/WorkingDirectory/OrthoFinder/Results_Oct17_2/WorkingDirectory/"
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

rule scatter_fasta:
    input:
        "resource/6_BLASTp/{prefix}.fasta"
    output:
        "resource/6_BLASTp/partition/{prefix}_{i}.fasta"
    params:
        n_partitions=n_partitions_blastp,
        i_partition=lambda w: int(w.i)
    conda:
        "env/hinflata.yaml"
    script:
        "scripts/scatter_fasta.py"

rule blastp:
    input:
         query="resource/6_BLASTp/partition/{prefix}.fasta",
         db="/data/zeynep/databases"
    output:
          "output/3_BLASTp/{prefix}.blastp"
    params:
            #format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sskingdoms stitle staxids",
            format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids",
            num_threads=1,
            evalue=1e-5,
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
