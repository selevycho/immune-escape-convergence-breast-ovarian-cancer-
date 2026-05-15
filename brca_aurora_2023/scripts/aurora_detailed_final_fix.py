import pandas as pd
import os

# 1. Абсолютные пути к твоим реальным данным
HOME = "/home/fr/fr_fr/fr_os136"
DATA_DIR = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/data")
RES_DIR = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results")

# Файлы данных
MUT_FILE = os.path.join(DATA_DIR, "data_mutations.txt")
Z_FILE = os.path.join(RES_DIR, "line1_aurora_merged_data.csv")

# Списки генов
HIGH_LIST = os.path.join(RES_DIR, "high_silencing_aurora_all.csv")
LOW_LIST = os.path.join(RES_DIR, "low_silencing_aurora_all.csv")

def create_report(target_genes_file, output_name):
    print(f"\n[+] Генерируем отчет: {output_name}")
    
    # Загружаем мутации (MAF/TXT формат обычно через таб)
    print("    - Читаем data_mutations.txt...")
    muts = pd.read_csv(MUT_FILE, sep='\t', low_memory=False)
    muts = muts[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification']]
    
    # Загружаем Z-scores
    print("    - Читаем Z-scores...")
    z_df = pd.read_csv(Z_FILE)
    
    # Находим колонку с ID образца в Z-файле
    s_col = 'Tumor_Sample_Barcode' if 'Tumor_Sample_Barcode' in z_df.columns else 'Sample_ID'
    
    # Объединяем
    merged = pd.merge(z_df, muts, left_on=['Hugo_Symbol', s_col], right_on=['Hugo_Symbol', 'Tumor_Sample_Barcode'])
    
    # Берем только гены из нашего списка
    genes = pd.read_csv(target_genes_file)['gene'].tolist()
    final_subset = merged[merged['Hugo_Symbol'].isin(genes)].copy()
    
    # Оставляем только ЗАГЛУШЕННЫЕ (Z < -1.0)
    silenced = final_subset[final_subset['Z_score'] < -1.0]

    if silenced.empty:
        print("    [!] Нет заглушенных мутаций для этого списка генов.")
        return

    # Создаем Pivot Table как на скриншоте
    report = pd.crosstab(silenced['Hugo_Symbol'], silenced['Variant_Classification'])
    
    # Колонки строго по ТЗ
    expected_cols = [
        'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
        'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 
        'Splice_Site', 'Translation_Start_Site'
    ]
    
    for c in expected_cols:
        if c not in report.columns:
            report[c] = 0
            
    report = report[expected_cols]
    report['Total_Silenced'] = report.sum(axis=1)
    report = report.sort_values(by='Total_Silenced', ascending=False).reset_index()

    # Сохранение
    save_path = os.path.join(RES_DIR, output_name)
    report.to_csv(save_path, index=False)
    print(f"    [OK] Сохранено в: {save_path}")
    print(report.head(10).to_string(index=False))

if __name__ == "__main__":
    if os.path.exists(MUT_FILE) and os.path.exists(Z_FILE):
        create_report(HIGH_LIST, "detailed_silenced_high_aurora.csv")
        create_report(LOW_LIST, "detailed_silenced_low_aurora.csv")
    else:
        print("[-] Ошибка: Проверь наличие data_mutations.txt и line1_aurora_merged_data.csv")
