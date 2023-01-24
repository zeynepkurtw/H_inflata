from Bio import SeqIO

try:
    multi_fasta = snakemake.input[0]
    partition_fasta = snakemake.output[0]
    n_partitions = snakemake.params.n_partitions
    i_partition = snakemake.params.i_partition
except NameError:
    # testing
    multi_fasta = "data/LGT_search/fasta_out/ss_trepo.fasta"
    partition_fasta = f"data/LGT_search/blast_in/partition/partition.fasta"

    n_partitions = 32
    i_partition = 0

print(locals())
with open(multi_fasta, "r") as fasta_in:
    with open(partition_fasta, "w") as fasta_out:
        i, j = 0, 0
        print(fasta_in)
        for i, seq in enumerate(SeqIO.parse(fasta_in, "fasta")):
            print(seq)
            if i % n_partitions == i_partition:
                SeqIO.write([seq], fasta_out, "fasta")
                j += 1
            i += 1

print(f"{j} out {i} sequences saved in {partition_fasta}.")
