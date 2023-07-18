import pandas as pd
import matplotlib.pyplot as plt

# assuming you have the tRNAscan-SE output in a tab-separated file called "output.txt"
df = pd.read_csv("output/5_tRNAscan/HIN.tRNAscan", sep="\t", header="infer", skiprows=1)

# plot the distribution of scores
plt.hist(df["Score"])
plt.xlabel("Score")
plt.ylabel("Count")
plt.title("Distribution of tRNA scores")
plt.show()

# plot the count of each tRNA type
df["tRNA #"].value_counts().plot(kind='bar')
plt.xlabel("tRNA Type")
plt.ylabel("Count")
plt.title("Count of each tRNA type")
plt.show()
