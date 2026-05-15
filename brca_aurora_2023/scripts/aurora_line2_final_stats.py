import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import os

# Пути
HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
MERGED_DATA = os.path.join(HOME, "results/line1_aurora_merged_data.csv")
EXP_FILE = os.path.join(HOME, "data/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
OUT_CSV = os.path.join(HOME, "results/hla_quartile_comparison_final.csv")

def main():
    print("[+] Запуск финального пересчета Line 2...")
    
    # 1. Загрузка и расчет Silencing Ratio
    df = pd.read_csv(MERGED_DATA)
    df['Is_Silenced'] = np.where(df['Z_score'] < -1.0, 1, 0)
    
    patient_stats = df.groupby('Tumor_Sample_Barcode').agg(
        Total_Muts=('Is_Silenced', 'count'),
        Silenced_Muts=('Is_Silenced', 'sum')
    ).reset_index()
    patient_stats['Ratio'] = patient_stats['Silenced_Muts'] / patient_stats['Total_Muts']

    # 2. Деление на квартили (строго 25/25)
    q75 = patient_stats['Ratio'].quantile(0.75)
    q25 = patient_stats['Ratio'].quantile(0.25)
    
    def get_group(r):
        if r >= q75: return 'High Silencing'
        if r <= q25: return 'Low Silencing'
        return 'Middle'
    
    patient_stats['Group'] = patient_stats['Ratio'].apply(get_group)
    extremes = patient_stats[patient_stats['Group'] != 'Middle'].copy()

    # 3. Загрузка экспрессии с нормализацией имен генов
    exp = pd.read_csv(EXP_FILE, sep='\t')
    # Мапим гены, игнорируя регистр и тире
    exp['gene_upper'] = exp['Hugo_Symbol'].str.upper().str.replace('-', '')
    targets = {'HLAA': 'HLA-A', 'HLAB': 'HLA-B', 'HLAC': 'HLA-C', 'B2M': 'B2M'}
    
    subset_exp = exp[exp['gene_upper'].isin(targets.keys())].copy()
    subset_exp['Hugo_Symbol'] = subset_exp['gene_upper'].map(targets)
    
    # Перевод в длинный формат
    melted = subset_exp.drop(columns=['gene_upper']).melt(
        id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Z_Exp'
    ).dropna()

    # 4. Мерджим и считаем
    final = pd.merge(extremes[['Tumor_Sample_Barcode', 'Group']], melted, on='Tumor_Sample_Barcode')
    final.to_csv(OUT_CSV, index=False)

    print(f"\n{'Ген':<8} | {'Mean Low':<12} | {'Mean High':<12} | {'P-value'}")
    print("-" * 55)
    
    for gene in ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']:
        g_data = final[final['Hugo_Symbol'] == gene]
        low = g_data[g_data['Group'] == 'Low Silencing']['Z_Exp']
        high = g_data[g_data['Group'] == 'High Silencing']['Z_Exp']
        
        if len(low) > 0 and len(high) > 0:
            p = ttest_ind(low, high, equal_var=False)[1]
            print(f"{gene:<8} | {low.mean():<12.3f} | {high.mean():<12.3f} | {p:.4f}")
        else:
            print(f"{gene:<8} | No Data      | No Data      | N/A")

if __name__ == "__main__":
    main()
