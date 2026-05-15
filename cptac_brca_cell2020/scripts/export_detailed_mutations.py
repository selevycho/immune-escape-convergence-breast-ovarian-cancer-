import pandas as pd
import numpy as np

def main():
    print("Генерация детализированных таблиц по типам мутаций...")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Загрузка и строгая фильтрация
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    valid_mutations = [
        'Missense_Mutation', 'Nonsense_Mutation', 
        'Frame_Shift_Del', 'Frame_Shift_Ins', 
        'In_Frame_Del', 'In_Frame_Ins',
        'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation'
    ]
    mut_filtered = mut[mut['Variant_Classification'].isin(valid_mutations)]
    
    # ВАЖНО: Теперь мы сохраняем колонку 'Variant_Classification' (тип мутации)
    unique_muts = mut_filtered.drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])

    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')
    exp_clean = exp_melted.dropna(subset=['Expression_FPKM'])

    # 2. Объединение ДНК и РНК
    neo_df = pd.merge(unique_muts[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification']], exp_clean, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    
    # 3. Разделяем все данные на два лагеря: Silenced и Expressed
    neo_df['Is_Silenced'] = np.where(neo_df['Expression_FPKM'] < 1.0, 1, 0)
    silenced_df = neo_df[neo_df['Is_Silenced'] == 1]
    expressed_df = neo_df[neo_df['Is_Silenced'] == 0]

    # 4. Функция для создания красивой сводной таблицы
    def create_pivot_table(df, total_col_name):
        if df.empty:
            return pd.DataFrame()
        
        # crosstab автоматически считает, сколько раз каждый тип мутации встретился в каждом гене
        pivot = pd.crosstab(df['Hugo_Symbol'], df['Variant_Classification'])
        
        # Добавляем колонку ИТОГО
        pivot[total_col_name] = pivot.sum(axis=1)
        
        # Сортируем по ИТОГО от большего к меньшему и сбрасываем индекс, чтобы Hugo_Symbol стал обычной колонкой
        pivot = pivot.sort_values(by=total_col_name, ascending=False).reset_index()
        
        # Заполняем пустоты нулями, если мутации такого типа не было
        pivot = pivot.fillna(0)
        return pivot

    # 5. Генерируем таблицы
    silenced_pivot = create_pivot_table(silenced_df, 'Total_Silenced')
    expressed_pivot = create_pivot_table(expressed_df, 'Total_Expressed')

    # 6. Сохраняем в CSV
    sil_file = 'results/detailed_silenced_mutations.csv'
    exp_file = 'results/detailed_expressed_mutations.csv'
    
    silenced_pivot.to_csv(sil_file, index=False)
    expressed_pivot.to_csv(exp_file, index=False)

    print(f"✅ Успешно! Таблицы с детализацией по типам мутаций сохранены:")
    print(f" - {sil_file}")
    print(f" - {exp_file}")

if __name__ == "__main__":
    main()
