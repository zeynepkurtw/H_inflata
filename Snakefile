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
         #expand("output/3_BLASTp/{file}_{i_partition}.blastp",file=["ss_trepo", "ss_hin", "og_hin_trepo"], i_partition=range(n_partitions_blastp)),
         #expand("output/7_tRNAscan/{sp}.tRNAscan", sp=["HIN", "muris", "wb", "spiro"]),
         #expand("output/7_tRNAscan/sensitive_search/{sp}.cov.tRNAscan",sp=[ "muris", "wb", "spiro"]),
         #expand("output/8_earlGrey/{sp}", sp=["HIN", "muris", "wb", "spiro"]),
         #"output/6_bedtools/HIN.intergenic.bed",
         #expand("output/6_bedtools/{sp}.intergenic.bed", sp=["HIN", "carpe","muris", "spiro", "kbiala" ]) #wb is missing
        expand("output/14_barRNAp/{sp}.rrna.gff", sp=["HIN", "muris", "wb", "spiro"]),

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
          "output/3_BLA STp/{prefix}.blastp"
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

rule tRNAscan:
    input:
         genome="resource/7_tRNAscan/{sp}.fasta"
    params: threads=8
    output:
          tRNA="output/7_tRNAscan/{sp}.tRNAscan",
          stats="output/7_tRNAscan/{sp}.stats"
    conda: "env/hinflata.yaml"
    script: "scripts/tRNAscan.py"

rule tRNAscan_cov:
    input:
         genome="resource/7_tRNAscan/{sp}.fasta"
    params: threads=32
    output:
          tRNA="output/7_tRNAscan/sensitive_search/{sp}.cov.tRNAscan",
          stats="output/7_tRNAscan/sensitive_search/{sp}.cov.stats"
    conda: "env/hinflata.yaml"
    script: "scripts/tRNAscan_cov.py"

rule earlGrey:
    input:
         genome="resource/10_earlGrey/{genome}.fasta"
    output:
          directory("output/8_earlGrey/{genome}")
    params:
        threads=8,
        species="{genome}"
    script:
          "scripts/earlGrey.py"

"""rule create_index:
    input:
         genome="resource/11_bedtools/HIN.fasta",
    output:
          index_genome="output/6_bedtools/HIN.fasta.fai",
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/create_index.py"


rule extract_intergenic:
    input:
         index_genome="output/6_bedtools/HIN.fasta.fai",
         genes="resource/11_bedtools/HIN_single-exon_79332.gff"
    output:
          cut_genome="output/6_bedtools/HIN_cut.fasta.fai",
          genes_sorted="output/6_bedtools/gff.sorted",
          intergenic="output/6_bedtools/HIN.intergenic.bed"
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/create_intergenic.py"
          """

rule create_index_:
    input:
         genome="resource/11_bedtools/{sp}.fasta",
    output:
          index_genome="output/6_bedtools/{sp}.fasta.fai",
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/create_index.py"


rule extract_intergenic_:
    input:
         index_genome="output/6_bedtools/{sp}.fasta.fai",
         genes="resource/11_bedtools/{sp}.CDS.gff"
    output:
          cut_genome="output/6_bedtools/{sp}.cut.fasta.fai",
          genes_sorted="output/6_bedtools/{sp}.CDS.sorted.gff",
          intergenic="output/6_bedtools/{sp}.intergenic.bed"
    conda:
         "env/hinflata.yaml"
    script:
          "scripts/create_intergenic.py"

rule barrnap:
    input:
        genome = "resource/14_barRNAp/{genome}.fasta"
    output:
        gff = "output/14_barRNAp/{genome}.rrna.gff",
        fasta = "resource/14_barRNAp/{genome}.rrna.fasta"

    conda: "env/hinflata.yaml"

    script:
          "scripts/barRNAp/barrnap.py"
