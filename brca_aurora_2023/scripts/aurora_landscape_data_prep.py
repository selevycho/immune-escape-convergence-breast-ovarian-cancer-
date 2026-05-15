import pandas as pd
import numpy as np
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
MUTS_FILE = os.path.join(HOME, "results/line1_aurora_merged_data.csv")
CNA_FILE = os.path.join(HOME, "data/data_log2_cna.txt")
OUTPUT_CSV = os.path.join(HOME, "results/aurora_patient_landscape_data.csv")

def main():
    print("[+] Шаг 1 & 2: Сборка таблицы по пациентам...")

    # 1. Обработка мутаций (Line 1)
    df_muts = pd.read_csv(MUTS_FILE)
    # Выделяем ID пациента: AUR_01_11_01 -> AUR_01_11
    df_muts['Patient_ID'] = df_muts['Tumor_Sample_Barcode'].str.split('_').str[:3].str.join('_')
    
    # Считаем Silencing Ratio
    df_muts['is_silenced'] = np.where(df_muts['Z_score'] < -1.0, 1, 0)
    # Группируем по пациенту
    pat_silencing = df_muts.groupby('Patient_ID').agg(
        ratio=('is_silenced', 'mean'),
        total_muts=('is_silenced', 'count')
    ).reset_index()

    # 2. Обработка CNA (Line 3)
    # Entrez IDs: B2M=567, HLA-A=3105, HLA-C=3107
    cna = pd.read_csv(CNA_FILE, sep='\t', low_memory=False)
    target_ids = [567, 3105, 3107]
    cna_subset = cna[cna['Entrez_Gene_Id'].isin(target_ids)].copy()
    
    # Маппинг для удобства
    id_map = {567: 'B2M', 3105: 'HLA_A', 3107: 'HLA_C'}
    cna_subset['gene'] = cna_subset['Entrez_Gene_Id'].map(id_map)
    
    # Melt (из широкого в длинный формат)
    cna_melted = cna_subset.drop(columns=['Entrez_Gene_Id']).melt(id_vars=['gene'], var_name='Sample_ID', value_name='log2_cna')
    cna_melted['Patient_ID'] = cna_melted['Sample_ID'].str.split('_').str[:3].str.join('_')
    
    # Группируем по пациенту и гену (берем MIN log2 - самый сильный делеционный сигнал)
    pat_cna = cna_melted.groupby(['Patient_ID', 'gene'])['log2_cna'].min().unstack().reset_index()
    
    # Создаем колонки статуса делеции (1 - удален, 0 - цел)
    for gene in ['B2M', 'HLA_A', 'HLA_C']:
        pat_cna[f'{gene}_status'] = np.where(pat_cna[gene] < -0.3, 1, 0)

    # 3. Финальное объединение (только те, у кого есть и РНК, и CNA)
    final_df = pd.merge(pat_silencing, pat_cna, on='Patient_ID', how='inner')
    
    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"[OK] Таблица на {len(final_df)} уникальных пациентов сохранена в {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
