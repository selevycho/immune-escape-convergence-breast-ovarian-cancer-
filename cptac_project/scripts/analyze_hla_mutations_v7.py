import pandas as pd
import os
import gzip
import glob

def main():
    # 1. Загружаем метаданные ПРАВИЛЬНО
    # Нам нужны колонки 'File ID' (имя папки) и 'Project ID'
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    
    # Создаем словарь маппинга: Folder_Name -> Project_Name
    # Важно: убираем возможные пробелы
    mapping = dict(zip(meta['File ID'].str.strip(), meta['Project ID'].str.strip()))
    
    # Создаем функцию для определения рака по Case ID (самый надежный способ)
    case_mapping = dict(zip(meta['File ID'].str.strip(), meta['Case ID'].str.strip()))

    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    maf_files = glob.glob(os.path.join('data/', '**/*.maf.gz'), recursive=True)
    
    all_results = []
    print(f"Total files to scan: {len(maf_files)}")

    for path in maf_files:
        # Извлекаем имя папки (UUID файла)
        folder_id = os.path.basename(os.path.dirname(path)).strip()
        
        # Пытаемся найти проект
        project = mapping.get(folder_id, "Unknown")
        case_id = case_mapping.get(folder_id, "Unknown")
        
        # Определяем рак по Case ID
        cancer_type = "Unknown"
        if 'BR' in str(case_id): cancer_type = 'Breast_Cancer'
        elif 'OV' in str(case_id): cancer_type = 'Ovarian_Cancer'
        elif 'CO' in str(case_id): cancer_type = 'Colon_Cancer'

        try:
            with gzip.open(path, 'rt') as f:
                header = None
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                if not header: continue
                
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    subset['Final_Project'] = cancer_type
                    all_results.append(subset)
        except:
            continue

    if not all_results:
        print("No mutations found.")
        return

    final = pd.concat(all_results)
    
    print("\n" + "="*60)
    print("ULTIMATE MUTATION REPORT (By Folder-to-Case mapping)")
    print("="*60)
    
    # Таблица 1
    print("\n[Table 1] All HLA/B2M Mutations by Real Cancer Type:")
    print(pd.crosstab(final['Final_Project'], final['Hugo_Symbol']))
    
    # Таблица 2
    lof = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = final[final['Variant_Classification'].isin(lof)]
    print("\n[Table 2] Severe LoF Mutations (The Killers):")
    print(pd.crosstab(dangerous['Final_Project'], dangerous['Hugo_Symbol']))
    
    final.to_csv('results/hla_mutations_v7.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
