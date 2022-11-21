import pandas as pd
from Bio import SeqIO
import glob
import seaborn as sns

path_interpro = "/Users/zeyku390/PycharmProjects/H_inflata/jupyter/data/interpro_annot/*.csv"
path_fasta = '/Users/zeyku390/PycharmProjects/data/proteome/*.fa'

""" Read files """
def read_files(path):
    list_files = glob.glob(path)
    dic = {}
    for element in list_files:
        i = element.split(".")[0]
        i = i.split("/")[-1]
        dic[i] = element
    return dic

""" Get gene ids from fasta """
def get_id(fasta_file):
    list_id = []
    with open(fasta_file, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            list_id.append(record.id)
    return  pd.DataFrame(list_id)

""" Substract functional genes from hypothetical ones """
def anti_join(x, y):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    # ans= ans.reset_index()
    return ans

"""def get_hypo(dic,sp):
    dic_hypo = {}
    dic_func = {}

    dic_func[sp] = dic[sp]
    dic_func[sp]["type"] = "functional"
    dic_hypo[sp] = anti_join(dic_id[sp], dic_func[sp])
    dic_hypo[sp]["type"] = "hypothetical"
    return dic_hypo, dic_func
"""


""" Gene ids from fasta """
dic_id = {}
for sp, fasta in read_files(path_fasta).items():
    dic_id[sp] = get_id(fasta)

""" Interpro annotation """
dic_ann = {}
for sp_ann, interpro in read_files(path_interpro).items():
    dic_ann[sp_ann] = pd.read_csv(interpro, sep="\t", header=None).drop(columns=[1,2]).drop_duplicates(0)


"""Functional and hypothetical proteins    """
dic_hypo = {}
dic_func = {}
dic_upset = {}
dic_count = {}
for sp_ann, ann in dic_ann.items():
    sp = sp_ann.split("_")[0]
    dic_func[sp] = dic_ann[sp_ann]
    dic_func[sp]["type"] = "functional"

    dic_hypo[sp] = anti_join(dic_id[sp], dic_func[sp])
    dic_hypo[sp]["type"] = "hypothetical"

    dic_upset[sp]=pd.concat([dic_hypo[sp], dic_func[sp]], axis=0)
    dic_upset[sp]["sp"] = sp
    dic_count[sp] = dic_upset[sp]["type"].value_counts().reset_index()

    df_upset = pd.concat(dic_upset, axis=0)
    df_count= pd.concat(dic_count, axis=0).reset_index().sort_values(by="type" ,ascending=False)


""" Bar plot    """
sns.set_style("whitegrid", {'axes.grid': False})
pal = sns.dark_palette("#69d", n_colors=2, reverse=False)
ax = sns.barplot(x='level_0', y='type', hue='index', data=df_count, palette=pal)
ax.legend().set_title("Interproscan annotations")
ax.set_xlabel('')
ax.set_ylabel('Number of genes')
ax.set_xticklabels(['H.inflata', 'K.bialata', 'C.membranifera', 'Trepomonas pc1', 'S.salmonicida', 'G.intestinalis', 'G.muris'], rotation=45, style='italic',ha="right")

plot = ax.get_figure()
plot.show()
plot.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1a.png', format="png",  bbox_inches='tight', dpi=1200)

