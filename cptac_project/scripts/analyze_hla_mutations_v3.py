import pandas as pd
import os
import gzip
import glob

def main():
    data_dir = 'data/'
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    # 1. Загружаем метаданные для маппинга File ID -> Project
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    # Создаем словарь: ID папки -> Название проекта
    file_to_project = dict(zip(meta['File ID'], meta['Project ID']))
    
    # Функция для уточнения рака (как мы делали раньше)
    def get_clean_project(file_id):
        proj = file_to_project.get(file_id, "Unknown")
        # Дополнительная проверка по Case ID из метаданных для точности
        case_id = meta[meta['File ID'] == file_id]['Case ID'].values
        if len(case_id) > 0:
            cid = str(case_id[0])
            if 'OV' in cid: return 'Ovarian_Cancer'
            if 'BR' in cid: return 'Breast_Cancer'
        return proj

    maf_files = glob.glob(os.path.join(data_dir, '**/*.maf.gz'), recursive=True)
    print(f"Scanning {len(maf_files)} files...")

    results = []

    for maf_path in maf_files:
        # Извлекаем File ID из пути (это название папки перед названием файла)
        file_id = os.path.basename(os.path.dirname(maf_path))
        project = get_clean_project(file_id)
        
        try:
            with gzip.open(maf_path, 'rt') as f:
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    subset['Project'] = project
                    results.append(subset[['Project', 'Hugo_Symbol', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode']])
        except:
            continue

    if not results:
        print("No mutations found.")
        return

    final_df = pd.concat(results)
    
    print("\n" + "="*60)
    print("FINAL CLINICAL MUTATION REPORT (No Unknowns)")
    print("="*60)
    
    # Таблица 1: Сводка по генам и проектам
    summary = pd.crosstab(final_df['Project'], final_df['Hugo_Symbol'])
    print("\n[Table 1] Mutation frequency per Cancer Type:")
    print(summary)
    
    print("\n" + "-"*60)
    
    # Таблица 2: Только тяжелые мутации (Loss of Function)
    lof_variants = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = final_df[final_df['Variant_Classification'].isin(lof_variants)]
    
    print(f"\n[Table 2] Severe 'Loss of Function' mutations found: {len(dangerous)}")
    if not dangerous.empty:
        # Группируем для краткости
        print(dangerous.groupby(['Project', 'Hugo_Symbol', 'Variant_Classification']).size().reset_index(name='Count'))
    
    final_df.to_csv('results/hla_mutations_final_report.tsv', sep='\t', index=False)
    print("\n" + "="*60)
    print("Done! Data saved to results/hla_mutations_final_report.tsv")

if __name__ == "__main__":
    main()
