#%%
import pandas as pd
import glob

" Extract  genes from eggnog"

#Read data
path = '/Users/zeyku390/PycharmProjects/comp_genomics/resources/eggnog/*annotations'
list_files = glob.glob(path)

#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element

#create dictinary for each species and its interproscan result
dic = {}
for key, value in sp_dic.items():
    dic[key] = ipr = pd.read_csv(value, header="infer", sep="\t", skiprows=3)

df_HIN_cn = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_HIN_cn.tsv",
                          sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})
df_HIN_ss = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_HIN_ss.tsv",
                          sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})

df_spiro_cn = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_spiro_cn.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})
df_spiro_ss = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_spiro_ss.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})

df_wb_cn = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_wb_cn.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})
df_wb_ss = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_wb_ss.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})

df_muris_cn = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_muris_cn.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})
df_muris_ss = pd.read_csv("/Users/zeyku390/PycharmProjects/comp_genomics/jupyter/ortofinder_4_sp/df_muris_ss.tsv",
                        sep="\t").drop_duplicates().rename(columns={"0": 0, "12": 12})
"HIN"

df_HIN_ss = pd.merge(df_HIN_ss, dic["HIN_NOG"], right_on="#query_name", left_on=0)[[0, "COG Functional cat."]]
df_HIN_ss = df_HIN_ss.drop_duplicates(subset=0)

cog_ss = df_HIN_ss.groupby("COG Functional cat.").count()
#cog_ss=cog_ss[cog_ss > 10].dropna()
cog_ss = cog_ss.reset_index().sort_values(by=0, ascending=False)
cog_ss

df_HIN_cn = pd.merge(df_HIN_cn, dic["HIN_NOG"], right_on="#query_name", left_on=0)[[0, "COG Functional cat."]]
df_HIN_cn = df_HIN_cn.drop_duplicates(subset=0)

cog_cn = df_HIN_cn.groupby("COG Functional cat.").count()
#cog_cn=cog_cn[cog_cn > 10].dropna()
cog_cn = cog_cn.reset_index().sort_values(by=0, ascending=False)
cog_cn

cog_ss = cog_ss.rename(columns={0: "species_specific"})
cog_ss
cog_cn = cog_cn.rename(columns={0: "conserved"})
cog_cat = pd.merge(cog_ss, cog_cn, how="outer").rename(columns={"COG Functional cat.": "COG cat"})
cog_cat

data_HIN = cog_cat.melt('COG cat', var_name="a", value_name='b' )
data_HIN["sp"]="HIN"
data_HIN =data_HIN.dropna(how="any")#.sort_values(by="b", ascending=False)
data_HIN
#%%
df_spiro_ss=pd.merge(df_spiro_ss, dic["spiro_NOG"], right_on="#query_name", left_on=0)[[0,"COG cat"]]
df_spiro_ss=df_spiro_ss.drop_duplicates(subset=0)

cog_ss=df_spiro_ss.groupby("COG cat").count()
#cog_ss=cog_ss[cog_ss > 100].dropna()
cog_ss=cog_ss.reset_index().sort_values(by=0, ascending=False)
cog_ss

df_spiro_cn=pd.merge(df_spiro_cn, dic["spiro_NOG"], right_on="#query_name", left_on=0)[[0, "COG cat"]]
df_spiro_cn=df_spiro_cn.drop_duplicates(subset=0)

cog_cn=df_spiro_cn.groupby("COG cat").count()
#cog_cn=cog_cn[cog_cn > 100].dropna()
cog_cn=cog_cn.reset_index().sort_values(by=0, ascending=False)
cog_cn

cog_ss=cog_ss.rename(columns={0:"species_specific"})
cog_ss
cog_cn=cog_cn.rename(columns={0:"conserved"})
cog_cat = pd.merge(cog_ss, cog_cn, how="outer")
cog_cat

data_spiro = cog_cat.melt('COG cat', var_name='a', value_name='b')
data_spiro["sp"]="spiro"
data_spiro =data_spiro.dropna(how="any")#.sort_values(by="b", ascending=False)
data_spiro
#%%
df_wb_ss=pd.merge(df_wb_ss, dic["wb_NOG"], right_on="#query_name", left_on=0)[[0,"COG cat"]]
df_wb_ss=df_wb_ss.drop_duplicates(subset=0)

cog_ss=df_wb_ss.groupby("COG cat").count()
#cog_ss=cog_ss[cog_ss > 100].dropna()
cog_ss=cog_ss.reset_index().sort_values(by=0, ascending=False)
cog_ss

df_wb_cn=pd.merge(df_wb_cn, dic["wb_NOG"], right_on="#query_name", left_on=0)[[0, "COG cat"]]
df_wb_cn=df_wb_cn.drop_duplicates(subset=0)

cog_cn=df_wb_cn.groupby("COG cat").count()
#cog_cn=cog_cn[cog_cn > 100].dropna()
cog_cn=cog_cn.reset_index().sort_values(by=0, ascending=False)
cog_cn

cog_ss=cog_ss.rename(columns={0:"species_specific"})
cog_ss
cog_cn=cog_cn.rename(columns={0:"conserved"})
cog_cat = pd.merge(cog_ss, cog_cn, how="outer")
cog_cat

data_wb = cog_cat.melt('COG cat', var_name='a', value_name='b')
data_wb["sp"]="wb"
data_wb =data_wb.dropna(how="any")#.sort_values(by="b", ascending=False)
data_wb
#%%
df_muris_ss=pd.merge(df_muris_ss, dic["muris_NOG"], right_on="#query_name", left_on=0)[[0,"COG cat"]]
df_muris_ss=df_muris_ss.drop_duplicates(subset=0)

cog_ss=df_muris_ss.groupby("COG cat").count()
#cog_ss=cog_ss[cog_ss > 100].dropna()
cog_ss=cog_ss.reset_index().sort_values(by=0, ascending=False)
cog_ss

df_muris_cn=pd.merge(df_muris_cn, dic["muris_NOG"], right_on="#query_name", left_on=0)[[0, "COG cat"]]
df_muris_cn=df_muris_cn.drop_duplicates(subset=0)

cog_cn=df_muris_cn.groupby("COG cat").count()
#cog_cn=cog_cn[cog_cn > 100].dropna()
cog_cn=cog_cn.reset_index().sort_values(by=0, ascending=False)
cog_cn

cog_ss=cog_ss.rename(columns={0:"species_specific"})
cog_ss
cog_cn=cog_cn.rename(columns={0:"conserved"})
cog_cat = pd.merge(cog_ss, cog_cn, how="outer")
cog_cat

data_muris = cog_cat.melt('COG cat', var_name='a', value_name='b')
data_muris["sp"]="muris"
data_muris =data_muris.dropna(how="any")#.sort_values(by="b", ascending=False)
data_muris
#%%
data_cat=pd.concat([data_HIN, data_spiro, data_wb, data_muris], axis=0)
data_cat
#%%
data_cat.groupby(["COG cat","a"]).first()
#%%
""""
A   RNA processing and modification
B	Chromatin Structure and dynamics
C	Energy production and conversion
D	Cell cycle control and mitosis
E	Amino Acid metabolis and transport
F	Nucleotide metabolism and transport
G	Carbohydrate metabolism and transport
H	Coenzyme metabolis
I	Lipid metabolism
J	Tranlsation
K	Transcription
L	Replication and repair
M	Cell wall/membrane/envelop biogenesis
N	Cell motility
O	Post-translational modification, protein turnover, chaperone functions
P	Inorganic ion transport and metabolism
Q	Secondary Structure
T	Signal Transduction
U	Intracellular trafficing and secretion
Y	Nuclear structure
Z	Cytoskeleton
R	General Functional Prediction only
S	Function Unknown
"""

import seaborn as sns
sns.set(rc={'figure.figsize':(30,5)})
sns.set_style("whitegrid", {'axes.grid' : False})
#pal=sns.dark_palette("#69d", n_colors=4)
#pal=sns.color_palette("ch:start=.2,rot=-.3", n_colors=4)
pal=sns.color_palette("Set2")

ax=sns.barplot(x='COG cat', y='b', hue='sp', data=data_cat.groupby("a").get_group("species_specific"), palette=pal)
#ax.legend( loc='upper right', borderaxespad=0).set_title("species specific COG")
ax.legend().set_title("species specific COG")
ax.set_ylabel('')
ax.set_xlabel('')
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

plot=ax.get_figure()
plot.savefig('/Users/zeyku390/PycharmProjects/comp_genomics/plots/barplot_COG_all_sp_ss.svg', format="svg",bbox_inches='tight')
#%%
import seaborn as sns
sns.set(rc={'figure.figsize':(30,5)})
sns.set_style("whitegrid", {'axes.grid' : False})
#pal=sns.dark_palette("#69d",n_colors=4)
pal=sns.color_palette("Set2")

ax=sns.barplot(x='COG cat', y='b', hue='sp', data=data_cat.groupby("a").get_group("conserved"), palette=pal)
ax.legend().set_title("Conserved COG")
ax.set_ylabel('')
ax.set_xlabel('')
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

plot=ax.get_figure()
plot.savefig('/Users/zeyku390/PycharmProjects/comp_genomics/plots/barplot_COG_all_sp_cn.svg', format="svg",bbox_inches='tight')
#%% md
# Subplots
#%%
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

#f,(ax1,ax2,ax3, ax4) = plt.subplots(4,1, gridspec_kw={'width_ratios':[1,1,1,1], "height_ratios":[1,1,1,1]} , figsize=(10,10))

f,(ax1,ax2,ax3, ax4) = plt.subplots(4,1, figsize=(15,15))

#f, (ax1,ax2,ax3, ax4) = plt.subplots(nrows=4, sharex=True, figsize=(30,10))
pal=sns.dark_palette("#69d", reverse=True ,n_colors=2)


ax_list = [ax1, ax2, ax3, ax4]
ax_list[0].get_shared_x_axes().join(ax_list[0], *ax_list[1:])

g1=sns.barplot(x='COG cat', y='b', hue='a', data=data_HIN[data_HIN["b"] >5], palette=pal, ax=ax1)
g2=sns.barplot(x='COG cat', y='b', hue='a', data=data_spiro[data_spiro["b"] >5], palette=pal, ax=ax2)
g3=sns.barplot(x='COG cat', y='b', hue='a', data=data_wb[data_wb["b"]> 5], palette=pal, ax=ax3)
g4=sns.barplot(x='COG cat', y='b', hue='a', data=data_muris[data_muris["b"] > 5], palette=pal, ax=ax4)



# may be needed to rotate the ticklabels correctly:
for ax in [g1,g2,g3,g4]:
    #tl = ax.get_xticklabels()
    #ax.set_xticklabels(tl, rotation=0)
    #tly = ax.get_yticklabels()
    #ax.set_yticklabels(tly, rotation=0)
    #ax.set_yticks([])
    #ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_label("")
    ax.legend([],[], frameon=False)


ax1.set_title("COG functional categories filtered genes>5")
ax1.legend(loc="upper right")

#plt=ax.get_figure()
plt.savefig('/Users/zeyku390/PycharmProjects/comp_genomics/plots/all_COG.svg', format='svg',bbox_inches='tight',dpi=1200)
plt.show()
#%% md
# sort by category
#%%
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

#f,(ax1,ax2,ax3, ax4) = plt.subplots(4,1, gridspec_kw={'width_ratios':[1,1,1,1], "height_ratios":[1,1,1,1]} , figsize=(10,10))

f,(ax1,ax2,ax3, ax4) = plt.subplots(4,1, figsize=(15,15))

#f, (ax1,ax2,ax3, ax4) = plt.subplots(nrows=4, sharex=True, figsize=(30,10))
pal=sns.dark_palette("#69d", reverse=True ,n_colors=2)


ax_list = [ax1, ax2, ax3, ax4]
ax_list[0].get_shared_x_axes().join(ax_list[0], *ax_list[1:])

"""order_list = ['first', 'second', 'third']
sns.barplot(x=df['x'], y=df['y'], order=order_list)"""

cog_order=["D", "M", "N", "O", "T", "U", "V", "W", "Y", "Z",
           "A", "B", "J", "K", "L",
           "C", "E", "F", "G", "H", "I", "P", "Q",
           "R", "S"]
g1=sns.barplot(x='COG cat', y='b', hue='a', data=data_HIN[data_HIN["b"] >5], color="darkseagreen", ax=ax1, order=cog_order)
g2=sns.barplot(x='COG cat', y='b', hue='a', data=data_spiro[data_spiro["b"] >5],color= "salmon", ax=ax2, order=cog_order)
g3=sns.barplot(x='COG cat', y='b', hue='a', data=data_wb[data_wb["b"]> 5], color= "salmon", ax=ax3, order=cog_order)
g4=sns.barplot(x='COG cat', y='b', hue='a', data=data_muris[data_muris["b"] > 5],color= "salmon", ax=ax4, order=cog_order)



# may be needed to rotate the ticklabels correctly:
for ax in [g1,g2,g3,g4]:
    #tl = ax.get_xticklabels()
    #ax.set_xticklabels(tl, rotation=0)
    #tly = ax.get_yticklabels()
    #ax.set_yticklabels(tly, rotation=0)
    #ax.set_yticks([])
    #ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_label("")
    ax.legend([],[], frameon=False)

ax1.set_ylabel("HIN")
ax2.set_ylabel("spiro")
ax3.set_ylabel("wb")
ax4.set_ylabel("muris")



ax1.set_title("COG functional categories filtered genes>5")
ax1.legend(loc="upper left")

#plt=ax.get_figure()
plt.savefig('/Users/zeyku390/PycharmProjects/H_inflata/plots/barplot_COG/all_COG_sorted_category.png', format='png',bbox_inches='tight',dpi=800)
plt.show()
