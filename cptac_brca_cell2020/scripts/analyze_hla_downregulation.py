import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def main():
    print("Запуск квартильного анализа HLA (Top 25% vs Bottom 25%)...")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Загрузка и фильтрация (Ultra-Robust)
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    valid_muts = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation']
    mut = mut[mut['Variant_Classification'].isin(valid_muts)]
    unique_muts = mut.drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])

    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])
    
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')
    exp_clean = exp_melted.dropna(subset=['Expression_FPKM'])

    # 2. Объединяем
    neo_df = pd.merge(unique_muts[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_clean, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    neo_df['Is_Silenced'] = np.where(neo_df['Expression_FPKM'] < 1.0, 1, 0)

    # 3. Считаем Ratio для пациентов
    patient_stats = neo_df.groupby('Tumor_Sample_Barcode').agg(
        Total_Mutations=('Is_Silenced', 'count'),
        Silenced_Count=('Is_Silenced', 'sum')
    ).reset_index()
    patient_stats['Silencing_Ratio'] = patient_stats['Silenced_Count'] / patient_stats['Total_Mutations']

    # 4. Выделяем экстремумы
    q75 = patient_stats['Silencing_Ratio'].quantile(0.75)
    q25 = patient_stats['Silencing_Ratio'].quantile(0.25)

    extremes = patient_stats[(patient_stats['Silencing_Ratio'] >= q75) | (patient_stats['Silencing_Ratio'] <= q25)].copy()
    
    # ВОЗВРАЩАЕМ КЛАССИЧЕСКИЕ НАЗВАНИЯ
    extremes['Group'] = np.where(extremes['Silencing_Ratio'] >= q75, 'High Silencing', 'Low Silencing')
    
    print("\nРазмер групп для анализа:")
    print(extremes['Group'].value_counts())

    # 5. Достаем экспрессию витрины
    target_genes = ['HLAA', 'HLAB', 'HLAC', 'B2M']
    hla_exp = exp_melted[exp_melted['Hugo_Symbol'].isin(target_genes)]
    
    final_df = pd.merge(extremes[['Tumor_Sample_Barcode', 'Group']], hla_exp, on='Tumor_Sample_Barcode', how='inner')
    final_df.to_csv('results/hla_comparison_data.csv', index=False)

    # 6. Статистика (t-test)
    print("\n=== СРАВНЕНИЕ ЭКСПРЕССИИ MHC-I (High vs Low Silencing) ===")
    print(f"{'Gene':<8} | {'Mean Low Silencing':<20} | {'Mean High Silencing':<20} | {'P-value':<10}")
    print("-" * 75)

    for gene in target_genes:
        gene_data = final_df[final_df['Hugo_Symbol'] == gene]
        low_grp = gene_data[gene_data['Group'] == 'Low Silencing']['Expression_FPKM'].dropna()
        high_grp = gene_data[gene_data['Group'] == 'High Silencing']['Expression_FPKM'].dropna()
        
        stat, pval = ttest_ind(low_grp, high_grp, equal_var=False)
        print(f"{gene:<8} | {low_grp.mean():<20.2f} | {high_grp.mean():<20.2f} | {pval:.4e}")

if __name__ == "__main__":
    main()
