import pandas as pd

def main():
    # Загружаем нашу итоговую таблицу метаданных
    df = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    
    # Создаем сводную таблицу: 
    # Строки - Тип ткани, Колонки - Тип данных (RNA-seq или Мутации)
    report = pd.crosstab(df['Tissue Type'], df['Data Type'])
    
    print("\n--- DATA AVAILABILITY REPORT ---")
    print(report)
    print("\n--------------------------------")
    
    # Ищем конкретно файлы экспрессии для нормы
    normal_rna = df[(df['Tissue Type'].str.contains('Normal', na=False)) & 
                    (df['Data Type'] == 'Gene Expression Quantification')]
    
    if not normal_rna.empty:
        print(f"BINGO! Found {len(normal_rna)} RNA-Seq files for Normal tissue.")
        print(normal_rna[['File ID', 'File Name', 'Case ID']].head())
    else:
        print("CONFIRMED: No RNA-Seq files found for Normal tissue in your downloads.")

if __name__ == "__main__":
    main()
