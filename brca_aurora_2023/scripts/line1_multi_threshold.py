import pandas as pd
import numpy as np

# 1. Загрузка данных
print("Loading AURORA 2023 data...")
mut_path = "../data/data_mutations.txt"
rna_path = "../data/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"

mutations = pd.read_csv(mut_path, sep='\t', low_memory=False, 
                        usecols=['Hugo_Symbol', 'Tumor_Sample_Barcode'])

rna = pd.read_csv(rna_path, sep='\t')
rna = rna.set_index('Hugo_Symbol')

# Список порогов для тестирования
thresholds = [-2.0, -1.5, -1.0, -0.5]

# 2. Подготовка к объединению
print("Mapping mutations to Z-scores...")
results_list = []

# Оптимизация: берем только те гены и образцы, которые есть в обоих датасетах
common_genes = list(set(mutations['Hugo_Symbol']) & set(rna.index))
rna_filtered = rna.loc[common_genes]

for idx, row in mutations.iterrows():
    gene = row['Hugo_Symbol']
    sample = row['Tumor_Sample_Barcode']
    
    if gene in rna_filtered.index and sample in rna.columns:
        z_val = rna.at[gene, sample]
        if not pd.isna(z_val):
            results_list.append({'gene': gene, 'z': z_val})

df_all = pd.DataFrame(results_list)
total_muts = len(df_all)

print(f"\nTotal Valid Matches Found: {total_muts}")
print("-" * 50)
print(f"{'Threshold (Z <)':<15} | {'Silenced Muts':<15} | {'Ratio (%)':<10}")
print("-" * 50)

# 3. Цикл по порогам
comparison_data = []
for t in thresholds:
    silenced_count = len(df_all[df_all['z'] < t])
    ratio = (silenced_count / total_muts) * 100
    print(f"{t:<15} | {silenced_count:<15} | {ratio:.2f}%")
    comparison_data.append({'Threshold': t, 'Ratio': ratio})

print("-" * 50)

# 4. Биологическая интерпретация (Топ гены для порога -1.0 как самого сбалансированного)
print("\nDeep Dive into Z < -1.0 (Realistic Low Expression):")
df_all['is_low'] = df_all['z'] < -1.0
gene_stats = df_all.groupby('gene')['is_low'].agg(['count', 'sum'])
gene_stats.columns = ['Total_Mut', 'Low_Exp_Mut']
gene_stats['Rate'] = (gene_stats['Low_Exp_Mut'] / gene_stats['Total_Mut']) * 100

print("\nTop 5 'Stealth' Genes at Z < -1.0:")
print(gene_stats[gene_stats['Total_Mut'] > 3].sort_values('Rate', ascending=False).head(5))

# 5. Сохранение итогов
pd.DataFrame(comparison_data).to_csv("../results/line1_multi_threshold_summary.csv", index=False)
print("\nMulti-threshold summary saved to results/ folder.")
