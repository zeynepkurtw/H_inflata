import pandas as pd
from utils import *


def group_stats(df):
    print("Number of OGs = ", len(df))
    for col in df.columns:
        if col != "Orthogroup":
            print(f'{col}\n')
            print(f"Sum of {col}: {df[col].sum()}")


def stats(df):
    df_sum= df.agg("sum", axis=0)
    return df_sum

df = get_og_count()

og_diplo = df[(df.loc[:, ["HIN", "trepo", "spiro", "wb", "muris"]] >= 1).all(1) &
              (df.loc[:, ["carpe", "kbiala"]] == 0).all(1)]

og_hin_trepo = df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                  (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)]

og_fl = df[(df.loc[:, ["carpe", "kbiala", "HIN", "trepo"]] >= 1).all(1) &
           (df.loc[:, ["spiro", "wb", "muris"]] == 0).all(1)]

og_hin_trepo_carpe = df[(df.loc[:, ["HIN", "trepo", "carpe"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "kbiala"]] == 0).all(1)]

og_hin_trepo_kbiala = df[(df.loc[:, ["HIN", "trepo", "kbiala"]] >= 1).all(1) &
                         (df.loc[:, ["spiro", "wb", "muris", "carpe"]] == 0).all(1)]


sum_diplo = stats(og_diplo)
sum_hin_trepo = stats(og_hin_trepo)
sum_fl = stats(og_fl)
sum_hin_trepo_carpe = stats(og_hin_trepo_carpe)
sum_hin_trepo_kbiala = stats(og_hin_trepo_kbiala)

#group 1 : og_diplo
gr1 = get_og_genes(og_diplo)



