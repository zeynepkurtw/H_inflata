import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

cont_blast= "resource/15_assembly_cont/prok_cont_hybrid_polished_masurca.blastn"
out_before= "plots/assembly_cont/Masurca_3_refprok_blastn_evalue_10.png"


data= pd.read_csv(cont_blast, header=None,sep='\t')

data=data.sort_values(by=[3], ascending=False)#.drop_duplicates(0)

df=data.iloc[:,[0,2,3,6,7,10,12,14]]
df=df.rename({0:"qseqid", 2:"pident", 3:"length", 6:"qstart",
              7:"qend", 10:"e_value" ,12:"qlen", 14:"stitle" },
            axis=1)
df1=df['stitle'].str.split(" ", expand=True)
df1['stitle_short'] = df1[[0,1]].agg(' '.join, axis=1)
df1=df1['stitle_short']
df=pd.concat([df,df1], axis=1).drop(columns="stitle")


import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
sns.set(rc={'figure.figsize':(5,5)})
sns.set_style("whitegrid", {'axes.grid' : False})

df1= df[df.e_value < pow(10,-10)]

plot2=sns.scatterplot(data=df1, x="qlen", y="length", hue="pident")
#place legend outside top right corner of plot
plot2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='pi')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel("contig")
plt.ylabel("alignment")
#plt.title( "Masurca_3refprok_blastn_evalue_minus10")
plt.savefig(out_before, dpi=300, bbox_inches='tight')
plt.show()

