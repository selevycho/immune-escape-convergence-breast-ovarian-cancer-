import pandas as pd
import os

def main():
    # Загружаем сырые данные мутаций
    if not os.path.exists('results/hla_mutations_raw.tsv'):
        print("Error: hla_mutations_raw.tsv not found. Run extract_hla_mutations.py first.")
        return
    
    muts = pd.read_csv('results/hla_mutations_raw.tsv', sep='\t')
    
    # Загружаем метаданные, чтобы понять где какой рак
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    
    # Уникальный маппинг Case ID к Проекту
    case_to_project = meta[['Case ID', 'Project ID']].drop_duplicates()
    
    # В MAF колонке Tumor_Sample_Barcode первые символы обычно соответствуют Case ID
    # Давай сделаем сопоставление по частичному совпадению
    def get_project(barcode):
        for _, row in case_to_project.iterrows():
            if str(row['Case ID']) in str(barcode):
                return row['Project ID']
        return "Unknown"

    print("Mapping mutations to projects...")
    muts['Project'] = muts['Tumor_Sample_Barcode'].apply(get_project)
    
    # Уточняем типы рака (Breast vs Ovarian) как мы делали раньше
    def refine_project(row):
        barcode = str(row['Tumor_Sample_Barcode'])
        if 'OV' in barcode: return 'Ovarian_Cancer'
        if 'BR' in barcode: return 'Breast_Cancer'
        return row['Project']

    muts['Project_Refined'] = muts.apply(refine_project, axis=1)

    print("\n" + "="*50)
    print("DETAILED MUTATION REPORT")
    print("="*50)
    
    # Сводка по проектам и генам
    summary = pd.crosstab(muts['Project_Refined'], muts['Hugo_Symbol'])
    print("\nMutation counts per gene and cancer type:")
    print(summary)
    
    print("\n" + "-"*50)
    print("TOP DANGEROUS MUTATIONS (Loss of Function):")
    # Ищем мутации, которые гарантированно ломают белок
    lof_variants = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = muts[muts['Variant_Classification'].isin(lof_variants)]
    
    if not dangerous.empty:
        print(dangerous[['Project_Refined', 'Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode']])
    else:
        print("No high-impact Loss-of-Function mutations found.")
    
    print("="*50)
    
    # Сохраняем детальный отчет
    muts.to_csv('results/hla_mutations_detailed.tsv', sep='\t', index=False)
    print(f"\nDetailed analysis saved to results/hla_mutations_detailed.tsv")

if __name__ == "__main__":
    main()
