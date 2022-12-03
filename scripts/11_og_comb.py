import pandas as pd
from utils import *

df = get_og_count()

""" OG Combinations """
og_comb = {
    "og_diplo": df[(df.loc[:, ["HIN", "trepo", "spiro", "wb", "muris"]] >= 1).all(1) &
                   (df.loc[:, ["carpe", "kbiala"]] == 0).all(1)],

    "og_hin_trepo" : df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)],

    "og_hin_trepo" : df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)],

    "og_fl" : df[(df.loc[:, ["carpe", "kbiala", "HIN", "trepo"]] >= 1).all(1) &
                 (df.loc[:, ["spiro", "wb", "muris"]] == 0).all(1)],

    "og_hin_trepo_kbiala" : df[(df.loc[:, ["HIN", "trepo", "kbiala"]] >= 1).all(1) &
                               (df.loc[:, ["spiro", "wb", "muris", "carpe"]] == 0).all(1)],

    "og_hin_trepo_carpe" : df[(df.loc[:, ["HIN", "trepo", "carpe"]] >= 1).all(1) &
                              (df.loc[:, ["spiro", "wb", "muris", "kbiala"]] == 0).all(1)]
}


def get_groups(og_comb, file_paths):
    dic_groups = {}

    for keys, values in og_comb.items():
        dic_groups[keys] = get_ann(get_og_genes(values))
        dic_groups[keys].to_csv(file_paths.format(keys), sep="\t", index=False)

    return dic_groups

def sum_stats():
    dic_stats = {}

    for keys, values in og_comb.items():
        dic_stats[keys] = values.agg("sum", axis=0)

    return dic_stats

groups = get_groups(og_comb, "data/og_groups/{}.csv")
group_stats = sum_stats()

