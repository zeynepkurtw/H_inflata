import pandas as pd
import matplotlib.pyplot as plt

#read ann_egg_cat file
og_ann=pd.read_csv("/Users/zeyku390/PycharmProjects/H_inflata/data/orthogroups/og_ann.csv", sep="\t", header='infer')

#drop nan values on COG column
og_ann_cog=og_ann[["COG","id", "sp" ]].dropna()
# groupby COG and sp
og_ann_cog=og_ann_cog.groupby(["COG"]).count().reset_index()


