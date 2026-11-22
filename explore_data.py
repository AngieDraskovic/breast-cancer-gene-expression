import pandas as pd
import matplotlib.pyplot as plt

file_path = 'data/c6a90de3-979b-46b0-9f84-2875c1e38742.rna_seq.augmented_star_gene_counts.tsv'

df = pd.read_csv(
    file_path,
    sep='\t',
    comment='#'
)

df = df[~df['gene_id'].str.startswith('N_')]

# print(df[['gene_name', 'unstranded']].head())
# print(df.shape) # (number of rows, number of cols)
# print(df.columns.tolist()) 

top_genes = df.nlargest(10, 'unstranded')

print(top_genes[['gene_name', 'unstranded']])

top_tpm = df.nlargest(10, 'tpm_unstranded')

print(top_tpm[['gene_name', 'tpm_unstranded']])



print(set(df['gene_type']))

print(df['gene_type'].value_counts())

print(df.groupby('gene_type').size().sort_values(ascending=False))

print(df['gene_type'].value_counts().plot(
    kind='bar',
    figsize=(14, 6)
))

# plt.show()
print(df.shape)
print(df.columns.tolist())