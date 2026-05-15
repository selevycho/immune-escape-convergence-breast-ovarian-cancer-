import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import os

# 1. Пути
HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
MERGED_DATA = os.path.join(HOME, "results/line1_aurora_merged_data.csv")
EXP_FILE = os.path.join(HOME, "data/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")

def main():
    print("[+] Запуск обновленного анализа HLA (фикс NaN)...")
    
    if not os.path.exists(MERGED_DATA):
        print(f"[-] Ошибка: Файл {MERGED_DATA} не найден.")
        return

    # 2. Загружаем данные по мутациям
    df = pd.read_csv(MERGED_DATA)
    df['Is_Silenced'] = np.where(df['Z_score'] < -1.0, 1, 0)

    # 3. Silencing Ratio
    patient_stats = df.groupby('Tumor_Sample_Barcode').agg(
        Total_Mutations=('Is_Silenced', 'count'),
        Silenced_Count=('Is_Silenced', 'sum')
    ).reset_index()
    patient_stats['Silencing_Ratio'] = patient_stats['Silenced_Count'] / patient_stats['Total_Mutations']

    # 4. Квартили
    q75 = patient_stats['Silencing_Ratio'].quantile(0.75)
    q25 = patient_stats['Silencing_Ratio'].quantile(0.25)
    
    def categorize(ratio):
        if ratio >= q75: return 'High Silencing'
        if ratio <= q25: return 'Low Silencing'
        return 'Middle'

    patient_stats['Group'] = patient_stats['Silencing_Ratio'].apply(categorize)
    extremes = patient_stats[patient_stats['Group'] != 'Middle'].copy()

    # 5. Загружаем экспрессию и фиксим имена генов
    exp = pd.read_csv(EXP_FILE, sep='\t')
    
    # Убираем дефисы и приводим к верхнему регистру для поиска
    exp['Search_Name'] = exp['Hugo_Symbol'].str.replace('-', '').str.replace(' ', '').str.upper()
    
    target_search = ['HLAA', 'HLAB', 'HLAC', 'B2M']
    hla_exp = exp[exp['Search_Name'].isin(target_search)].copy()
    
    # Возвращаем нормальные имена для красоты
    name_map = {'HLAA': 'HLA-A', 'HLAB': 'HLA-B', 'HLAC': 'HLA-C', 'B2M': 'B2M'}
    hla_exp['Hugo_Symbol'] = hla_exp['Search_Name'].map(name_map)

    # Melt
    hla_melted = hla_exp.drop(columns=['Search_Name']).melt(
        id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Z_score_HLA'
    )
    
    # Чистим NaN из экспрессии ДО объединения
    hla_melted = hla_melted.dropna(subset=['Z_score_HLA'])

    # 6. Объединение
    final_df = pd.merge(extremes[['Tumor_Sample_Barcode', 'Group']], hla_melted, on='Tumor_Sample_Barcode')
    final_df.to_csv(os.path.join(HOME, 'results/hla_quartile_comparison_aurora.csv'), index=False)

    # 7. Статистика
    print("\n=== СРАВНЕНИЕ MHC-I (Top 25% vs Bottom 25%) ===")
    print(f"{'Gene':<8} | {'Mean Low':<15} | {'Mean High':<15} | {'P-value'}")
    print("-" * 65)

    for gene in ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']:
        gene_data = final_df[final_df['Hugo_Symbol'] == gene]
        low_grp = gene_data[gene_data['Group'] == 'Low Silencing']['Z_score_HLA']
        high_grp = gene_data[gene_data['Group'] == 'High Silencing']['Z_score_HLA']
        
        if len(low_grp) > 1 and len(high_grp) > 1:
            stat, pval = ttest_ind(low_grp, high_grp, equal_var=False)
            print(f"{gene:<8} | {low_grp.mean():<15.3f} | {high_grp.mean():<15.3f} | {pval:.4e}")
        else:
            print(f"{gene:<8} | {'No Data':<15} | {'No Data':<15} | {'N/A'}")

if __name__ == "__main__":
    main()
