import pandas as pd
import numpy as np

def main():
    print("Запуск геномного анализа (Линия Обороны 3: Физическое удаление)...\n")
    cna_file = 'data/data_cna.txt'

    # 1. Загрузка данных CNA
    print(f"Чтение файла: {cna_file}")
    cna = pd.read_csv(cna_file, sep='\t')
    
    if 'Entrez_Gene_Id' in cna.columns:
        cna = cna.drop(columns=['Entrez_Gene_Id'])

    # 2. Фильтрация только нужных нам генов витрины
    target_genes = ['HLAA', 'HLAB', 'HLAC', 'B2M']
    cna_hla = cna[cna['Hugo_Symbol'].isin(target_genes)]

    if cna_hla.empty:
        print("Внимание: точные имена HLA не найдены, проверь формат файла...")
        return

    # 3. Переводим таблицу в удобный формат
    cna_melted = cna_hla.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='CNA_Value')
    cna_melted = cna_melted.dropna(subset=['CNA_Value'])
    
    # 4. Классификация стратегии выживания
    # Значения <= -1.0 означают потерю гетерозиготности (LOH) или полное удаление
    cna_melted['Is_Deleted'] = np.where(cna_melted['CNA_Value'] <= -1.0, 1, 0)

    # 5. Подсчет статистики
    stats = cna_melted.groupby('Hugo_Symbol').agg(
        Total_Patients=('Tumor_Sample_Barcode', 'count'),
        Deleted_Count=('Is_Deleted', 'sum')
    ).reset_index()

    stats['Deletion_Percent'] = (stats['Deleted_Count'] / stats['Total_Patients']) * 100
    stats = stats.sort_values(by='Deletion_Percent', ascending=False)

    # 6. Вывод результатов
    print("=== ЧАСТОТА ФИЗИЧЕСКОГО УДАЛЕНИЯ ВИТРИНЫ (Genomic Deletion) ===")
    print(f"{'Gene':<8} | {'Total Patients':<15} | {'Deleted (LOH/-2)':<18} | {'Percent':<10}")
    print("-" * 60)
    
    for index, row in stats.iterrows():
        print(f"{row['Hugo_Symbol']:<8} | {int(row['Total_Patients']):<15} | {int(row['Deleted_Count']):<18} | {row['Deletion_Percent']:.1f}%")

    # Сохраняем результат
    stats.to_csv('results/hla_genomic_deletions.csv', index=False)
    print("\n✅ Данные сохранены в results/hla_genomic_deletions.csv")

if __name__ == "__main__":
    main()
