import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('results/hla_matrix_v3.tsv', sep='\t')
m = pd.melt(df, id_vars=['Project'], value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'])

plt.figure(figsize=(12, 7))
# Явно задаем цвета, чтобы точно видеть разницу
palette = {'Breast_Cancer': 'seagreen', 'Ovarian_Cancer': 'orange', 'Colon_Cancer': 'dodgerblue', 'Other': 'gray'}

sns.boxplot(data=m, x='variable', y='value', hue='Project', palette=palette)
plt.yscale('log')
plt.title('HLA Expression by Cancer Type (CPTAC-2)')
plt.ylabel('TPM (Log Scale)')
plt.savefig('results/hla_comparison_v3.png')
print("Plot v3 saved.")
