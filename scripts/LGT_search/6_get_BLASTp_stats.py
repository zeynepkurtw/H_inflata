import pandas as pd
import glob


try:
    blast_file = snakemake.input.blast_file
    out = snakemake.output[0]
except NameError:
    # testing
    blast_file = "data/LGT_search/blast_out/*_filtered.blastp"

    out_file = "data/LGT_search/blast_out/{}_stats.blastp"


def kingdom_stats(group):

    grouped_df = group.groupby("")

    # Calculate the sum of the number of occurrences
    occurrence_sum = grouped_df['counts'].sum()

    # Calculate the percentage of occurrence per group
    grouped_df['percentage'] = grouped_df['counts'] / occurrence_sum * 100

    # Return the result as a dataframe
    return grouped_df

blast_files = glob.glob(blast_file)
stats = {}

for element in blast_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    stats[i] = element

for key, values in stats.items():

    df = pd.read_csv(values, sep='\t', header= None)

    grouped = df.groupby([0])

    stats[key] = grouped[16].apply(lambda x: x.value_counts(normalize=True)).round(2)

for key, values in stats.items():
    values = values.reset_index().rename(columns={ 0:"id", "level_1":"taxa", 16:"perc_taxa"})
    values.to_csv(out_file.format(key), sep="\t", index= False)