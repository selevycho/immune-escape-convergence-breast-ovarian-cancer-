import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load data
df = pd.read_csv('results/hla_matrix_v2.tsv', sep='\t')

# 2. Reshape for plotting
m = pd.melt(df, id_vars=['Project'], value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
            var_name='Gene', value_name='TPM')

# 3. Create a professional plot
plt.figure(figsize=(12, 7))
sns.set_style("whitegrid")

# Create boxplot with dots (stripplot) on top
sns.boxplot(data=m, x='Gene', y='TPM', hue='Project', palette='Set2', showfliers=False)
sns.stripplot(data=m, x='Gene', y='TPM', hue='Project', dodge=True, alpha=0.3, jitter=True, color='black', size=2)

plt.yscale('log') # Log scale is essential for gene expression
plt.title('HLA Class I Expression: Breast (BRCA) vs Ovarian (OV) Cancer', fontsize=15)
plt.ylabel('Expression (TPM) - Log Scale', fontsize=12)
plt.xlabel('Gene', fontsize=12)

# Save
plt.savefig('results/hla_comparison_v2.png', dpi=300)
print("SUCCESS: New plot saved to results/hla_comparison_v2.png")
