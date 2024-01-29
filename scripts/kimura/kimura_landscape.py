import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# read RepeatMasker .out file
# skip the first three rows, use whitespace as a separator
repeat_df = pd.read_csv('resource/12_kimura_landscape/HIN.fasta.out', sep='\s+', skiprows=3, header=None, error_bad_lines=False)

# extract the repeat class (5th column) and Kimura divergence (2nd column)
repeat_df = repeat_df[[4, 1]]

# calculate frequency of each repeat class-Kimura value pair
repeat_freq = repeat_df.groupby([4, 1]).size().reset_index()

# pivot the DataFrame to get repeat classes as columns and Kimura values as index
repeat_freq_pivot = repeat_freq.pivot(index=1, columns=4)

# calculate cumulative sums for each repeat class column
repeat_freq_cumsum = repeat_freq_pivot.cumsum()

# plot the repeat landscape
plt.stackplot(repeat_freq_cumsum.index, repeat_freq_cumsum.values.T)

plt.xlabel('Kimura divergence')
plt.ylabel('Cumulative repeat length (bp)')
plt.title('Repeat Landscape')
plt.show()
