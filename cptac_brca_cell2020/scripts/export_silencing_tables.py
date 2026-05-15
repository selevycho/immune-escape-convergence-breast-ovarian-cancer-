import pandas as pd
import numpy as np

def main():
    print("Генерация Топ-20 таблиц (High vs Low Silencing)...")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Загрузка и фильтрация
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    valid_mutations = [
        'Missense_Mutation', 'Nonsense_Mutation', 
        'Frame_Shift_Del', 'Frame_Shift_Ins', 
        'In_Frame_Del', 'In_Frame_Ins',
        'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation'
    ]
    mut_filtered = mut[mut['Variant_Classification'].isin(valid_mutations)]
    unique_muts = mut_filtered.drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])

    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')
    exp_clean = exp_melted.dropna(subset=['Expression_FPKM'])

    # 2. Объединение и расчет статуса
    neo_df = pd.merge(unique_muts[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_clean, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    neo_df['Is_Silenced'] = np.where(neo_df['Expression_FPKM'] < 1.0, 1, 0)

    # 3. Группировка
    gene_stats = neo_df.groupby('Hugo_Symbol').agg(
        Total_Mutations=('Tumor_Sample_Barcode', 'count'),
        Silenced_Count=('Is_Silenced', 'sum')
    ).reset_index()

    gene_stats['Expressed_Count'] = gene_stats['Total_Mutations'] - gene_stats['Silenced_Count']
    gene_stats['Silenced_Percent'] = (gene_stats['Silenced_Count'] / gene_stats['Total_Mutations']) * 100

    # 4. Фильтр шума (минимум 3 мутации)
    gene_stats_filtered = gene_stats[gene_stats['Total_Mutations'] >= 3]

    # 5. Выделяем Топ-20 High Silencing (Невидимки)
    top_high = gene_stats_filtered.sort_values(by=['Silenced_Percent', 'Total_Mutations'], ascending=[False, False]).head(20)
    
    # 6. Выделяем Топ-20 Low Silencing (Драйверы)
    top_low = gene_stats_filtered.sort_values(by=['Silenced_Percent', 'Total_Mutations'], ascending=[True, False]).head(20)

    # 7. Сохранение в файлы
    high_file = 'results/high_silencing_top20.csv'
    low_file = 'results/low_silencing_top20.csv'
    
    top_high.to_csv(high_file, index=False)
    top_low.to_csv(low_file, index=False)

    print(f"✅ Успешно! Таблицы сохранены в папку results:\n - {high_file}\n - {low_file}")

if __name__ == "__main__":
    main()
