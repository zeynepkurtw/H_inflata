import pandas as pd
import glob


"OG diplo subset"
df=pd.read_csv( "/Users/zeyku390/PycharmProjects/H.inflata/output/1_orthofinder/Results_Oct10_7sp/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", header='infer')
df=df[["Orthogroup","HIN", "trepo", "spiro", "wb", "muris", "kbiala", "carpe", "Total"]]

og_diplo=df[(df.iloc[:,[1,2,3,4,5]] >=1).all(1) & (df.iloc[:,[6,7]] ==0).all(1)]

"OG gene list"
OG_gene_list = pd.read_csv('/Users/zeyku390/PycharmProjects/H.inflata/output/1_orthofinder/Results_Oct10_7sp/Orthogroups/Orthogroups.txt', header=None, dtype=str, delim_whitespace=True)
OG_gene_list = OG_gene_list.apply(lambda x: x.str.replace(":", ""))


"OG diplo subset gene list with annotations"
def OG_stack(OG_list):
    df = pd.merge(OG_gene_list, OG_list["Orthogroup"], right_on="Orthogroup", left_on=0).drop(
        columns="Orthogroup").rename(columns={0:"OG"}).set_index("OG")
    df=df.stack().reset_index().drop(columns=["level_1"])
    return df

df_og_stack=OG_stack(og_diplo)

"Add annotations"
path = '/Users/zeyku390/PycharmProjects/H.inflata/jupyter/data/*.csv'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split("_annot")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element

dic_annot= {}
for key, value in sp_dic.items():
    dic_annot[key] = pd.read_csv(value, sep="\t", header="infer")

def add_annot(df, df_annot):
    df=pd.merge(df, df_annot, right_on="id", left_on=0, how="left").drop(columns=[0]).dropna()
    return df

df1=add_annot(OG_stack(og_diplo), dic_annot["HIN"])
df2=add_annot(OG_stack(og_diplo), dic_annot["spiro"])
df3=add_annot(OG_stack(og_diplo), dic_annot["muris"])
df4=add_annot(OG_stack(og_diplo), dic_annot["wb"])
df5=add_annot(OG_stack(og_diplo), dic_annot["trepo"])

df_og_stack_annot=pd.concat([df1,df2,df3,df4,df5], axis=0)
df_og_stack_annot.to_csv("/Users/zeyku390/PycharmProjects/H.inflata/jupyter/1_orthofinder/data/og_diplo_annot.csv",sep="\t", header="infer", index=False)
