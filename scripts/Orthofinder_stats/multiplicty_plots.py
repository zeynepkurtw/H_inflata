import pandas as pd
import numpy as np
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt

#input
path="/Users/zeyku390/PycharmProjects/comp_genomics/resources/3_orthofinder_4_sp/"
df1=pd.read_csv(path+"OrthologuesStats_one-to-one.tsv",
               sep="\t", header='infer',index_col="Unnamed: 0").rename(columns={"HIN_79341_annotation":"HIN"}, index={"HIN_79341_annotation":"HIN"})
df2=pd.read_csv(path+"OrthologuesStats_one-to-many.tsv",
               sep="\t", header='infer',index_col="Unnamed: 0").rename(columns={"HIN_79341_annotation":"HIN"}, index={"HIN_79341_annotation":"HIN"})
df3=pd.read_csv(path+"OrthologuesStats_many-to-one.tsv",
               sep="\t", header='infer',index_col="Unnamed: 0").rename(columns={"HIN_79341_annotation":"HIN"}, index={"HIN_79341_annotation":"HIN"})
df4=pd.read_csv(path+"OrthologuesStats_many-to-many.tsv",
               sep="\t", header='infer',index_col="Unnamed: 0").rename(columns={"HIN_79341_annotation":"HIN"}, index={"HIN_79341_annotation":"HIN"})

#output
path_out="/Users/zeyku390/PycharmProjects/H_inflata/plots/Orthofinder_4sp/"

df_hin=pd.concat([df1.iloc[0],df2.iloc[0],df3.iloc[0],df4.iloc[0]], axis=1)
df_hin.columns = ["1:1", "1:m", "m:1", "m:m"]

df_spiro=pd.concat([df1.iloc[2],df2.iloc[2],df3.iloc[2],df4.iloc[2]], axis=1)
df_spiro.columns = ["1:1", "1:m", "m:1", "m:m"]

df_muris=pd.concat([df1.iloc[1],df2.iloc[1],df3.iloc[1],df4.iloc[1]], axis=1)
df_muris.columns = ["1:1", "1:m", "m:1", "m:m"]

df_wb=pd.concat([df1.iloc[3],df2.iloc[3],df3.iloc[3],df4.iloc[3]], axis=1)
df_wb.columns = ["1:1", "1:m", "m:1", "m:m"]

df_sum=df1+ df2.values + df3.values + df4.values


def multiplicity_plot_hin(diplomonad,summation):
    df_sum_X=summation.iloc[0,:] #manipulate according to the organism
    df_sum_X = df_sum_X.to_numpy() #convert to numpy array
    df_sum_X = df_sum_X[:, np.newaxis]  # adds a new axis -> 2D array (6,) => (6,1)
    X=(diplomonad/df_sum_X)*100
    X=X.dropna()#.set_index(X.columns[0])

    sns.set(rc={'figure.figsize': (10, 5)})
    sns.set_style("whitegrid", {'axes.grid': False})
    pal = sns.diverging_palette( 150, 20, s=60, l=60, n=4, center="light")
    sns.set_palette(palette=pal)

    plot=X.plot.barh(stacked=True)
    #plt.title("H. inflata")
    plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.savefig(path_out + '/multiplicty_hin.png', format='png',
                bbox_inches='tight', dpi=800)
    plt.show()
    return X, plot


multiplicity_plot_hin(df_hin, df_sum)


def multiplicity_plot_spiro(diplomonad,summation):
    sns.set(rc={'figure.figsize': (10, 5)})
    sns.set_style("whitegrid", {'axes.grid': False})
    pal = sns.diverging_palette(150, 20, s=60, l=60, n=4, center="light")
    sns.set_palette(palette=pal)

    df_sum_X=summation.iloc[2,:] #manipulate according to the organism
    df_sum_X = df_sum_X.to_numpy() #convert to numpy array
    df_sum_X = df_sum_X[:, np.newaxis]  # adds a new axis -> 2D array (6,) => (6,1)
    X=(diplomonad/df_sum_X)*100
    X=X.dropna()#.set_index(X.columns[0])
    plot=X.plot.barh(stacked=True).legend(loc="best") #exclude the NaN values
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.savefig(path_out + '/multiplicty_spiro.png', bbox_inches='tight', dpi=800)
    plt.show()
    return X, plot

multiplicity_plot_spiro(df_spiro, df_sum)

def multiplicity_plot(diplomonad,summation):
    df_sum_X=summation.iloc[3,:] #manipulate according to the organism
    df_sum_X = df_sum_X.to_numpy() #convert to numpy array
    df_sum_X = df_sum_X[:, np.newaxis]  # adds a new axis -> 2D array (6,) => (6,1)
    X=(diplomonad/df_sum_X)*100
    X=X.dropna()#.set_index(X.columns[0])

    plot = X.plot.barh(stacked=True).legend(loc="best")  # exclude the NaN values
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.savefig(path_out + '/multiplicty_wb.png', bbox_inches='tight', dpi=800)
    plt.show()
    return X, plot

multiplicity_plot(df_wb, df_sum)

def multiplicity_plot(diplomonad,summation):
    df_sum_X=summation.iloc[1,:] #manipulate according to the organism
    df_sum_X = df_sum_X.to_numpy() #convert to numpy array
    df_sum_X = df_sum_X[:, np.newaxis]  # adds a new axis -> 2D array (6,) => (6,1)
    X=(diplomonad/df_sum_X)*100
    X=X.dropna()#.set_index(X.columns[0])

    plot = X.plot.barh(stacked=True).legend(loc="best")  # exclude the NaN values
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.savefig(path_out + '/multiplicty_muris.png', bbox_inches='tight', dpi=800)
    plt.show()
    return X, plot

multiplicity_plot(df_muris, df_sum)