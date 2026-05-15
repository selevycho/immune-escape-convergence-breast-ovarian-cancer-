import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import os

# 1. Пути
HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
Z_FILE = os.path.join(HOME, "data/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
MUT_FILE = os.path.join(HOME, "data/data_mutations.txt")
HIGH_GENES_FILE = os.path.join(HOME, "results/detailed_silenced_high_aurora.csv")
LOW_GENES_FILE = os.path.join(HOME, "results/detailed_silenced_low_aurora.csv")

def run_line2_analysis():
    print("[+] Loading AURORA datasets...")
    # Загружаем Z-scores экспрессии
    z_scores = pd.read_csv(Z_FILE, sep='\t')
    # Загружаем мутации
    mutations = pd.read_csv(MUT_FILE, sep='\t', low_memory=False)
    
    # Списки генов из Line 1
    high_genes = pd.read_csv(HIGH_GENES_FILE)['Hugo_Symbol'].unique()
    low_genes = pd.read_csv(LOW_GENES_FILE)['Hugo_Symbol'].unique()

    # 2. Определяем группы пациентов
    # Группа 1: Пациенты с мутациями в генах High Silencing
    patients_high = mutations[mutations['Hugo_Symbol'].isin(high_genes)]['Tumor_Sample_Barcode'].unique()
    # Группа 2: Пациенты с мутациями в генах Low Silencing
    patients_low = mutations[mutations['Hugo_Symbol'].isin(low_genes)]['Tumor_Sample_Barcode'].unique()

    # Убираем пересечения (если пациент в обеих группах, для чистоты анализа)
    # Но в метастазах часто есть и то, и то. Оставим как есть или возьмем уникальных.
    
    # 3. Извлекаем данные HLA/B2M
    hla_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    z_hla = z_scores[z_scores['Hugo_Symbol'].isin(hla_genes)].set_index('Hugo_Symbol').T
    z_hla.index.name = 'Tumor_Sample_Barcode'

    # Создаем итоговый датасет для сравнения
    z_hla = z_hla.dropna()
    
    def get_group(barcode):
        if barcode in patients_low: return 'Low_Silencing_Group'
        if barcode in patients_high: return 'High_Silencing_Group'
        return 'Other'

    z_hla['Group'] = z_hla.index.map(get_group)
    comparison_df = z_hla[z_hla['Group'] != 'Other'].copy()

    # 4. Статистика
    results = []
    print("\n" + "="*50)
    print(f"{'Gene':<10} | {'High Mean':<10} | {'Low Mean':<10} | {'P-Value'}")
    print("-" * 50)

    for gene in hla_genes:
        group_high = comparison_df[comparison_df['Group'] == 'High_Silencing_Group'][gene]
        group_low = comparison_df[comparison_df['Group'] == 'Low_Silencing_Group'][gene]
        
        stat, p = mannwhitneyu(group_high, group_low)
        
        results.append({
            'Gene': gene,
            'High_Mean': group_high.mean(),
            'Low_Mean': group_low.mean(),
            'P_value': p
        })
        
        print(f"{gene:<10} | {group_high.mean():<10.3f} | {group_low.mean():<10.3f} | {p:.4f}")

    # Сохраняем данные для построения графика
    comparison_df.to_csv(os.path.join(HOME, "results/line2_comparison_data.csv"))
    print("="*50)
    print(f"[+] Statistics saved and comparison data ready for plotting.")

if __name__ == "__main__":
    run_line2_analysis()
