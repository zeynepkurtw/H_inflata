import pandas as pd

df =pd.read_csv("/Users/zeyku390/Downloads/wikipedia_300.csv")

df['Target'] = df['Category'].apply(lambda c: 0 if 'Programming' in c else 1)
