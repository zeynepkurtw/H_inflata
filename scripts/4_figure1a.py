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
    annot= []
    with open(fasta_file, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            list_id.append(record.id)
            annot.append(record.description)
    df = pd.DataFrame(list(zip(list_id, annot)), columns=[0, "desc"])
    df["desc"] = df.apply(lambda x: x["desc"].replace(x[0], "").strip(), axis=1)
    return df #pd.DataFrame(list_id)

""" Substract functional genes from hypothetical ones """
def anti_join(x, y):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    # ans= ans.reset_index()
    return ans


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
dic_cons_hypo = {}
dic_func = {}
dic_bar = {}
dic_count = {}
for sp_ann, ann in dic_ann.items():
    sp = sp_ann.split("_")[0]
    dic_func[sp] = dic_ann[sp_ann]
    dic_func[sp]["type"] = "functional"

    dic_hypo[sp] = anti_join(dic_id[sp], dic_func[sp])
    dic_hypo[sp]["type"] = "hypothetical"

    #dic_cons_hypo[sp] = dic_id[sp][dic_id[sp]["desc"].str.contains("Conserved hypothetical")]
    #dic_cons_hypo[sp]["type"] = "conserved_hypothetical"

    #dic_bar[sp]=pd.concat([dic_hypo[sp], dic_func[sp],dic_cons_hypo[sp]], axis=0)
    dic_bar[sp] = pd.concat([dic_hypo[sp], dic_func[sp]], axis=0)
    dic_bar[sp]["sp"] = sp
    dic_count[sp] = dic_bar[sp]["type"].value_counts().reset_index()

    df_bar = pd.concat(dic_bar, axis=0)
    df_count= pd.concat(dic_count, axis=0).reset_index().sort_values(by="type" ,ascending=False)


""" Bar plot    """
sns.set_style("whitegrid", {'axes.grid': False})
#pal = sns.dark_palette("#69d", n_colors=2, reverse=False)

colors = [ "#2A4A5C",  "#7A93A2", "#41738F"]
pal = sns.set_palette(sns.color_palette(colors))

#, fill=False, hatch='/'
order = ["carpe","kbiala","HIN","trepo", "spiro","wb", "muris"]
#hue_order= ["functional", "hypothetical", "conserved_hypothetical"] # remove conserved hypothetical
hue_order= ["functional", "hypothetical"]
#ax1 = sns.barplot(x='level_0', y='type', hue='index', data=df_count, palette=pal, order=order, fill=False, edgecolor="#2A4A5C")
ax1 = sns.barplot(x='level_0', y='type', hue='index', data=df_count, palette=pal, order=order, hue_order=hue_order)


ax1.legend(loc='upper right').set_title("annotations" )
ax1.set_xlabel('')
ax1.set_ylabel('Number of genes')
ax1.set_xticklabels(['C. membranifera',  'K. bialata', 'H. inflata', 'Trepomonas pc1', 'S. salmonicida', 'G. intestinalis', 'G. muris'], rotation=45, style='italic',ha="right")

plot = ax1.get_figure()
plot.show()
plot.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/figure1a_.png', format="png",  bbox_inches='tight', dpi=1200)

