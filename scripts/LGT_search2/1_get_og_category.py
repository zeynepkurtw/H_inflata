import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann = "data/orthogroups/og_ann.csv"
    og_count_s = "data/orthogroups/og_count_s.csv"

    out_file_csv = "data/LGT_search/{}.csv"
    out_file_xlsx= "data/LGT_search/{}.xlsx"


def get_og_genes(og_ann, og_category: pd.DataFrame) -> pd.DataFrame:
    """
    @param og_category: Orthgroups with a species combination.
    @return: Orthogroup with its genes for given category.
    """
    og_ann_ = pd.read_csv(og_ann, sep="\t", header="infer")
    df = pd.merge(og_ann_, og_category)

    return df


df = pd.read_csv(og_count_s, sep="\t", header='infer')


""" OG Combinations """
og_comb = {

    "ss_hin" : df[(df.loc[:, ["HIN"]] >= 1).all(1) &
                            (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala", "trepo"]] == 0).all(1)],
    "ss_trepo": df[(df.loc[:, ["trepo"]] >= 1).all(1) &
                 (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala", "HIN"]] == 0).all(1)],

    "og_hin_trepo" : df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)],

}


dic_groups = {}
for keys, values in og_comb.items():

    dic_groups[keys] = get_og_genes(og_ann, values).sort_values(by=[ "OG", "Total"])

    dic_groups[keys].to_csv(out_file_csv.format(keys), sep="\t", index=False)
    dic_groups[keys].to_excel(out_file_xlsx.format(keys), index=False)


dic_stats={}
""" og_comb statistics """
for keys, values in og_comb.items():
    dic_stats[keys] = values.agg("sum", axis=0)

