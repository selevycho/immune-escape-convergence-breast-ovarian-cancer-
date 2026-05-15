import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import os

# Пути
HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
CNA_FILE = os.path.join(HOME, "data/data_log2_cna.txt")
GROUPS_FILE = os.path.join(HOME, "results/hla_quartile_comparison_final.csv")

# Маппинг Entrez ID -> Hugo Symbol
GENE_MAP = {
    567: 'B2M',
    3105: 'HLA-A',
    3106: 'HLA-B',
    3107: 'HLA-C'
}

def main():
    print("[+] Запуск Line 3: Анализ по Entrez ID (B2M и HLA)...")
    
    # 1. Загружаем группы
    groups_df = pd.read_csv(GROUPS_FILE).drop_duplicates(subset=['Tumor_Sample_Barcode'])
    
    # 2. Читаем CNA (только нужные строки для экономии памяти)
    print("    - Читаю data_log2_cna.txt...")
    # Читаем файл, так как он большой, отфильтруем Entrez_Gene_Id сразу
    cna = pd.read_csv(CNA_FILE, sep='\t', low_memory=False)
    
    # Оставляем только наши гены
    cna_targets = cna[cna['Entrez_Gene_Id'].isin(GENE_MAP.keys())].copy()
    cna_targets['Hugo_Symbol'] = cna_targets['Entrez_Gene_Id'].map(GENE_MAP)

    # 3. Транспонируем (переводим образцы в строки)
    cna_melted = cna_targets.drop(columns=['Entrez_Gene_Id']).melt(
        id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Log2_Ratio'
    )
    
    # 4. Гармонизация ID пациентов
    # В CNA: AUR_03_549-2_03 -> приводим к стандарту (замена _ на - и т.д. если нужно)
    # Но судя по head, они и так в AUR формате. Просто мерджим.
    final_df = pd.merge(groups_df[['Tumor_Sample_Barcode', 'Group']], cna_melted, on='Tumor_Sample_Barcode')
    
    final_df['Log2_Ratio'] = pd.to_numeric(final_df['Log2_Ratio'], errors='coerce')
    final_df = final_df.dropna(subset=['Log2_Ratio'])

    # Порог делеции (Log2 < -0.3)
    final_df['Is_Deletion'] = np.where(final_df['Log2_Ratio'] < -0.3, 1, 0)

    print("\n=== РЕЗУЛЬТАТЫ LINE 3 (AURORA) ===")
    print(f"{'Ген':<8} | {'Low Silencing (%)':<20} | {'High Silencing (%)':<20} | {'P-value'}")
    print("-" * 80)

    summary = []
    for gene in GENE_MAP.values():
        g_data = final_df[final_df['Hugo_Symbol'] == gene]
        if g_data.empty: continue
        
        counts = pd.crosstab(g_data['Group'], g_data['Is_Deletion'])
        # Проверка наличия колонок
        for col in [0, 1]:
            if col not in counts.columns: counts[col] = 0
            
        l_del, l_total = counts.loc['Low Silencing', 1], counts.loc['Low Silencing'].sum()
        h_del, h_total = counts.loc['High Silencing', 1], counts.loc['High Silencing'].sum()
        
        _, p = fisher_exact([[l_del, l_total-l_del], [h_del, h_total-h_del]])
        
        print(f"{gene:<8} | {l_del/l_total*100:<19.1f}% | {h_del/h_total*100:<19.1f}% | {p:.4f}")
        summary.append({'Gene': gene, 'Group': 'Low Silencing', 'Del_Freq': l_del/l_total})
        summary.append({'Gene': gene, 'Group': 'High Silencing', 'Del_Freq': h_del/h_total})

    out_path = os.path.join(HOME, "results/line3_deletion_summary.csv")
    pd.DataFrame(summary).to_csv(out_path, index=False)
    print(f"\n[OK] Данные сохранены в: {out_path}")

if __name__ == "__main__":
    main()
