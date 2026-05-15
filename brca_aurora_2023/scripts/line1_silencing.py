import pandas as pd
import numpy as np

# 1. Загрузка данных
print("Loading AURORA 2023 data...")
mut_path = "../data/data_mutations.txt"
rna_path = "../data/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"

# Загружаем мутации (берем только нужные колонки для экономии памяти)
mutations = pd.read_csv(mut_path, sep='\t', low_memory=False, 
                        usecols=['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification'])

# Загружаем RNA Z-scores (ID гена в первой колонке)
rna = pd.read_csv(rna_path, sep='\t')
rna = rna.set_index('Hugo_Symbol')

# 2. Очистка ID образцов
# В cBioPortal ID в РНК и Мутациях часто немного отличаются (например, -01 в конце).
# Приводим их к единому виду для сопоставления.
mutations['Sample_ID'] = mutations['Tumor_Sample_Barcode']
rna_samples = rna.columns[1:] # Пропускаем колонку Entrez_Gene_Id

# 3. Анализ Silencing (Линия 1)
print("Analyzing Transcriptomic Silencing (Z-score threshold < -2.0)...")

silenced_count = 0
total_valid_mutations = 0
gene_stats = []

# Проходим по всем мутациям и ищем их экспрессию
for idx, row in mutations.iterrows():
    gene = row['Hugo_Symbol']
    sample = row['Sample_ID']
    
    if gene in rna.index and sample in rna.columns:
        z_score = rna.at[gene, sample]
        
        # Проверяем на NaN (бывает, если данных по гену нет)
        if not pd.isna(z_score):
            total_valid_mutations += 1
            is_silenced = z_score < -2.0 # Твой новый порог для Z-scores
            
            if is_silenced:
                silenced_count += 1
            
            gene_stats.append({
                'Hugo_Symbol': gene,
                'Sample_ID': sample,
                'Z_score': z_score,
                'Is_Silenced': is_silenced
            })

# 4. Расчет результатов
stats_df = pd.DataFrame(gene_stats)
global_ratio = (silenced_count / total_valid_mutations) * 100

print("-" * 30)
print(f"GLOBAL RESULTS FOR AURORA 2023:")
print(f"Total mutations analyzed: {total_valid_mutations}")
print(f"Silenced mutations (Z < -2.0): {silenced_count}")
print(f"Mutation Silencing Ratio: {global_ratio:.2f}%")
print("-" * 30)

# 5. Анализ по генам (High vs Low Silencing)
gene_summary = stats_df.groupby('Hugo_Symbol').agg(
    Total_Mutations=('Is_Silenced', 'count'),
    Silenced_Mutations=('Is_Silenced', 'sum')
)
gene_summary['Silencing_Rate'] = (gene_summary['Silenced_Mutations'] / gene_summary['Total_Mutations']) * 100

# Вывод топ-10 самых "скрытных" генов (Passenger candidates)
print("\nTop 10 High-Silencing Genes (Stealth Targets):")
print(gene_summary[gene_summary['Total_Mutations'] > 2].sort_values('Silencing_Rate', ascending=False).head(10))

# Вывод топ-10 самых "открытых" генов (Driver candidates)
print("\nTop 10 Low-Silencing Genes (Exposed Drivers):")
print(gene_summary[gene_summary['Total_Mutations'] > 5].sort_values('Silencing_Rate', ascending=True).head(10))

# 6. Сохранение результатов для Jupyter
stats_df.to_csv("../results/line1_aurora_results.csv", index=False)
gene_summary.to_csv("../results/line1_aurora_gene_summary.csv")
print("\nResults saved to results/ folder.")
