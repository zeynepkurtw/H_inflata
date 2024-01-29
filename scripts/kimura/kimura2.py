import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

# Open the file
with open('resource/12_kimura_landscape/HIN.fasta.out', 'r') as f:
    lines = f.readlines()

# Remove the '*' from each line
lines = [line.rstrip(' *\n') for line in lines]

# Join the lines back into a single string
data = '\n'.join(lines)

# Use StringIO to read the string as a file
df = pd.read_csv(StringIO(data), sep='\s+', skiprows=3, header=None)

# Handle the overlap
grouped = df.groupby([10])[1].sum().reset_index()
grouped.columns = ['Kimura_Divergence', 'Count']

# Plot
plt.figure(figsize=(10,6))
plt.bar(grouped['Kimura_Divergence'], grouped['Count'], color='skyblue', edgecolor='blue')
plt.xlabel('Kimura Divergence')
plt.ylabel('Count')
plt.title('Kimura Divergence Landscape')
# Rotate the x labels by 45 degrees
plt.xticks(rotation=45)
# Adjust bottom margin
plt.subplots_adjust(bottom=0.15)

plt.show()