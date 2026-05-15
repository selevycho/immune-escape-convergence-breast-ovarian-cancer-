import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def main():
    print("Loading data...")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Загрузка мутаций
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    valid_mutations = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']
    mut = mut[mut['Variant_Classification'].isin(valid_mutations)][['Hugo_Symbol', 'Tumor_Sample_Barcode']]

    # 2. Загрузка экспрессии
    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])

    # --- ИСПРАВЛЕННЫЕ НАЗВАНИЯ ГЕНОВ (без дефисов!) ---
    target_genes = ['HLAA', 'HLAB', 'HLAC', 'B2M']
    actual_targets = [g for g in target_genes if g in exp['Hugo_Symbol'].unique()]
    print(f"Найдены гены для анализа: {actual_targets}\n")

    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')

    # 3. Считаем статус мутаций
    neo_df = pd.merge(mut, exp_melted, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    neo_df['Silenced'] = np.where(neo_df['Expression_FPKM'] < 1.0, 1, 0)

    # 4. Делим пациентов на 2 группы
    patient_stats = neo_df.groupby('Tumor_Sample_Barcode').agg(
        Total_Mutations=('Hugo_Symbol', 'count'),
        Silenced_Mutations=('Silenced', 'sum')
    ).reset_index()
    patient_stats['Silencing_Ratio'] = patient_stats['Silenced_Mutations'] / patient_stats['Total_Mutations']

    median_ratio = patient_stats['Silencing_Ratio'].median()
    patient_stats['Group'] = np.where(patient_stats['Silencing_Ratio'] >= median_ratio, 'High Silencing', 'Low Silencing')

    # 5. Достаем экспрессию и сохраняем
    hla_exp = exp_melted[exp_melted['Hugo_Symbol'].isin(actual_targets)]
    final_df = pd.merge(patient_stats[['Tumor_Sample_Barcode', 'Group']], hla_exp, on='Tumor_Sample_Barcode', how='inner')
    final_df.to_csv('results/hla_groups_expression.tsv', sep='\t', index=False)

    # 6. DEA
    print("--- ДИФФЕРЕНЦИАЛЬНАЯ ЭКСПРЕССИЯ (DEA) ---")
    print(f"{'Gene':<8} | {'Mean Low Silencing':<20} | {'Mean High Silencing':<20} | {'P-value':<10}")
    print("-" * 65)

    for gene in actual_targets:
        gene_data = final_df[final_df['Hugo_Symbol'] == gene]
        low_group = gene_data[gene_data['Group'] == 'Low Silencing']['Expression_FPKM'].dropna()
        high_group = gene_data[gene_data['Group'] == 'High Silencing']['Expression_FPKM'].dropna()
        
        if len(low_group) > 0 and len(high_group) > 0:
            stat, pval = ttest_ind(low_group, high_group, equal_var=False)
            print(f"{gene:<8} | {low_group.mean():<20.2f} | {high_group.mean():<20.2f} | {pval:.4e}")
        else:
            print(f"{gene:<8} | Недостаточно данных")

if __name__ == "__main__":
    main()
