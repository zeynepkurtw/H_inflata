import pandas as pd
import glob

"OG diplo subset"
df = pd.read_csv("output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv",
    sep="\t", header='infer')
df = df[["Orthogroup", "HIN", "trepo", "spiro", "wb", "muris", "kbiala", "carpe", "Total"]]

og_diplo = df[(df.iloc[:, [1, 2, 3, 4, 5]] >= 1).all(1) & (df.iloc[:, [6, 7]] == 0).all(1)]

"OG gene list"
og_gene_list = pd.read_csv(
    'output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt',
    header=None, dtype=str, delim_whitespace=True)
og_gene_list = og_gene_list.apply(lambda x: x.str.replace(":", ""))

"OG diplo subset gene list with annotations"


def og_stack(og_list):
    df = pd.merge(og_gene_list, og_list["Orthogroup"], right_on="Orthogroup", left_on=0).drop(
        columns="Orthogroup").rename(columns={0: "OG"}).set_index("OG")
    df = df.stack().reset_index().drop(columns=["level_1"])
    return df


df_og_stack = og_stack(og_diplo)

"Add annotations"
path = 'jupyter/data/*.csv'
list_files = glob.glob(path)

# get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split("_annot")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element

dic_annot = {}
for key, value in sp_dic.items():
    dic_annot[key] = pd.read_csv(value, sep="\t", header="infer")


def add_annot(df, df_annot):
    df = pd.merge(df, df_annot, right_on="id", left_on=0, how="left").drop(columns=[0]).dropna()
    return df


for value in dic_annot.values():
    df_og_stack_annot = pd.concat([add_annot(og_stack(og_diplo), value)], axis=0)

#df = pd.concat([add_annot(og_stack(og_diplo), value) for value in dic_annot.values()], axis=0) # fancy

df_og_stack_annot.to_csv("jupyter/1_orthofinder/data/og_diplo_annot.csv",
                         sep="\t", header=True, index=False)
