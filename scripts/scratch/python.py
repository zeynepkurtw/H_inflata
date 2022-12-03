from utils import *

def get_groups():
    #group 1 : og_diplo
    gr1 = get_og_genes(og_diplo)
    gr1 = get_ann(gr1)

    #group 2 : og_hin_trepo
    gr2 = get_og_genes(og_hin_trepo)
    gr2 = get_ann(gr2)

    # group 3 : og_fl
    gr3 = get_og_genes(og_fl)
    gr3 = get_ann(gr3)

    # group 4 : og_hin_trepo_kbiala
    gr4 = get_og_genes(og_hin_trepo_kbiala)
    gr4 = get_ann(gr4)

    # group 5 : og_hin_trepo_carpe
    gr5 = get_og_genes(og_hin_trepo_carpe)
    gr5 = get_ann(gr5)

    return gr1, gr2, gr3, gr4, gr5

def get_groups1():
    return [get_ann(get_og_genes(og)) for og in [og_diplo, og_hin_trepo, og_fl, og_hin_trepo_kbiala, og_hin_trepo_carpe]]

test1 = get_groups()

test2 = get_groups1()

test3 = []
for i in range(10):
    test3.append(2*i)

test4 = [2 * i for i in range(10)]


def group_stats(df):
    print("Number of OGs = ", len(df))
    for col in df.columns:
        if col != "Orthogroup":
            print(f'{col}\n')
            print(f"Sum of {col}: {df[col].sum()}")

def get_groups(og_comb, file_paths):
    dfs = []
    for og in og_comb:
        df = get_ann(get_og_genes(og))
        df.to_csv(file_paths.format(og), sep="\t", index=False)
        dfs.append(df)

    return dfs

def sum_stats():
    return [stats(keys) for keys, values in og_comb]