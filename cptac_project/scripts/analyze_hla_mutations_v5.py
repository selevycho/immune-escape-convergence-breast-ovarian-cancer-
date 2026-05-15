import pandas as pd
import os
import gzip
import glob

def main():
    # 1. Загружаем метаданные
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    
    # Делаем маппинг: File ID -> Case ID
    file_to_case = dict(zip(meta['File ID'], meta['Case ID']))
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    maf_files = glob.glob(os.path.join('data/', '**/*.maf.gz'), recursive=True)
    
    all_results = []
    print(f"Mapping {len(maf_files)} files using Folder IDs...")

    for path in maf_files:
        # ID папки - это имя родительской директории файла
        folder_id = os.path.basename(os.path.dirname(path))
        case_id = file_to_case.get(folder_id)
        
        if not case_id:
            continue
            
        # Определяем тип рака по Case ID
        cancer = "Unknown"
        if 'BR' in str(case_id): cancer = 'Breast_Cancer'
        elif 'OV' in str(case_id): cancer = 'Ovarian_Cancer'
        
        try:
            with gzip.open(path, 'rt') as f:
                # Пропускаем комменты до заголовка
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                if not subset.empty:
                    subset['Real_Cancer'] = cancer
                    subset['Case_ID'] = case_id
                    all_results.append(subset)
        except:
            continue

    if not all_results:
        print("Still nothing. Check if File IDs in metadata match data/ folder names.")
        return

    final = pd.concat(all_results)
    
    print("\n" + "="*60)
    print("FINAL ACCURATE REPORT (By Folder Mapping)")
    print("="*60)
    
    print("\nTotal Mutations per Gene and Cancer:")
    print(pd.crosstab(final['Real_Cancer'], final['Hugo_Symbol']))
    
    print("\n" + "-"*60)
    print("Severe (Loss of Function) Mutations:")
    lof = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = final[final['Variant_Classification'].isin(lof)]
    print(pd.crosstab(dangerous['Real_Cancer'], dangerous['Hugo_Symbol']))
    
    final.to_csv('results/hla_mutations_final_v5.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
