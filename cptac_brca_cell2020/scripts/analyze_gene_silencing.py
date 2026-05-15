import pandas as pd
import numpy as np

def main():
    print("Analyzing Gene-Specific Silencing Rates...\n")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Загрузка и очистка данных (тот же Ultra-Robust подход)
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

    # 2. Объединяем ДНК и РНК
    neo_df = pd.merge(unique_muts[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_clean, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')

    # 3. Определяем статус (1 = Заглушен, 0 = Экспрессируется)
    neo_df['Is_Silenced'] = np.where(neo_df['Expression_FPKM'] < 1.0, 1, 0)

    # 4. Группируем по имени гена
    gene_stats = neo_df.groupby('Hugo_Symbol').agg(
        Total_Mutations=('Tumor_Sample_Barcode', 'count'),
        Silenced_Count=('Is_Silenced', 'sum')
    ).reset_index()

    # Считаем производные метрики
    gene_stats['Expressed_Count'] = gene_stats['Total_Mutations'] - gene_stats['Silenced_Count']
    gene_stats['Silenced_Percent'] = (gene_stats['Silenced_Count'] / gene_stats['Total_Mutations']) * 100

    # 5. ФИЛЬТР: Оставляем только те гены, которые мутировали минимум у 3 пациентов
    gene_stats_filtered = gene_stats[gene_stats['Total_Mutations'] >= 3]

    # Сортируем: сначала те, что глушатся чаще всего (100%), затем по общему числу мутаций
    top_silenced = gene_stats_filtered.sort_values(by=['Silenced_Percent', 'Total_Mutations'], ascending=[False, False])
    
    # Сортируем наоборот: те, что НЕ глушатся (0% Silenced), это драйверы рака
    top_expressed = gene_stats_filtered.sort_values(by=['Silenced_Percent', 'Total_Mutations'], ascending=[True, False])

    # 6. Вывод в терминал
    print("=== ТОП-10 ГЕНОВ, КОТОРЫЕ РАК ПРЯЧЕТ ЧАЩЕ ВСЕГО (High Silencing) ===")
    print(top_silenced[['Hugo_Symbol', 'Total_Mutations', 'Silenced_Count', 'Silenced_Percent']].head(10).to_string(index=False))

    print("\n=== ТОП-10 ГЕНОВ, КОТОРЫЕ РАК НИКОГДА НЕ ПРЯЧЕТ (Low Silencing / Drivers) ===")
    print(top_expressed[['Hugo_Symbol', 'Total_Mutations', 'Expressed_Count', 'Silenced_Percent']].head(10).to_string(index=False))

    # 7. Сохраняем полную таблицу для отчета
    output_file = 'results/gene_specific_silencing.csv'
    top_silenced.to_csv(output_file, index=False)
    print(f"\nFull gene statistics saved to: {output_file}")

if __name__ == "__main__":
    main()
