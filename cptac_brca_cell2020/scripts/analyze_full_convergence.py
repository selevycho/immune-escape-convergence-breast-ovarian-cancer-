import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def main():
    print("📊 Запуск анализа CONVERGENCE на ПОЛНОЙ выборке (N=122)...")
    
    # 1. Загрузка данных для расчета Silencing Ratio (из 1-го пункта)
    mut = pd.read_csv('data/data_mutations.txt', sep='\t', low_memory=False)
    exp = pd.read_csv('data/data_mrna_seq_fpkm.txt', sep='\t')
    
    # Очистка мутаций
    valid_muts = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site']
    mut = mut[mut['Variant_Classification'].isin(valid_muts)].drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])
    
    # Очистка экспрессии
    if 'Entrez_Gene_Id' in exp.columns: exp = exp.drop(columns=['Entrez_Gene_Id'])
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')

    # Расчет Silencing Ratio для ВСЕХ
    merged = pd.merge(mut[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_melted, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    merged['Is_Silenced'] = np.where(merged['Expression_FPKM'] < 1.0, 1, 0)
    
    patient_stats = merged.groupby('Tumor_Sample_Barcode').agg(
        Total_Mutations=('Is_Silenced', 'count'),
        Silenced_Count=('Is_Silenced', 'sum')
    ).reset_index()
    patient_stats['Silencing_Ratio'] = patient_stats['Silenced_Count'] / patient_stats['Total_Mutations']

    # 2. Разделение на две группы по МЕДИАНЕ (50/50)
    median_val = patient_stats['Silencing_Ratio'].median()
    patient_stats['Group'] = np.where(patient_stats['Silencing_Ratio'] >= median_val, 'High Silencing', 'Low Silencing')

    # 3. Загрузка данных CNA (ДНК) по B2M
    cna = pd.read_csv('data/data_cna.txt', sep='\t')
    b2m_cna = cna[cna['Hugo_Symbol'] == 'B2M'].melt(
        id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='CNA_Value'
    )
    b2m_cna['B2M_Status'] = np.where(b2m_cna['CNA_Value'] <= -1, 'Deleted', 'Intact')

    # 4. Объединение
    final_df = pd.merge(patient_stats[['Tumor_Sample_Barcode', 'Group']], b2m_cna[['Tumor_Sample_Barcode', 'B2M_Status']], on='Tumor_Sample_Barcode')

    # 5. Анализ
    table = pd.crosstab(final_df['Group'], final_df['B2M_Status'])
    odds_ratio, p_value = fisher_exact(table)

    print(f"\n--- ИТОГОВАЯ ТАБЛИЦА (N={len(final_df)}) ---")
    print(table)
    print("-" * 35)

    for group in table.index:
        total = table.loc[group].sum()
        deleted = table.loc[group, 'Deleted']
        print(f"{group}: {deleted} из {total} пациентов ({ (deleted/total)*100 :.1f}%) имеют делецию B2M")

    print(f"\nP-value: {p_value:.4f}")
    
    # Сохранение
    final_df.to_csv('results/full_convergence_results.csv', index=False)
    print("\n✅ Полные результаты сохранены в results/full_convergence_results.csv")

if __name__ == "__main__":
    main()
