import pandas as pd
import os
import gzip
import glob

def main():
    data_dir = 'data/'
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    # 1. Загружаем метаданные для маппинга
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    
    # 2. Создаем "железный" маппинг Case_ID -> Реальный Рак
    def identify_cancer(case_id):
        cid = str(case_id).upper()
        if 'OV' in cid: return 'Ovarian_Cancer'
        if 'BR' in cid: return 'Breast_Cancer'
        return 'Unknown'

    # Создаем словарь: Case ID -> Тип рака
    case_to_cancer = {row['Case ID']: identify_cancer(row['Case ID']) for _, row in meta.iterrows()}

    maf_files = glob.glob(os.path.join(data_dir, '**/*.maf.gz'), recursive=True)
    results = []

    print(f"Scanning {len(maf_files)} files and mapping to correct cancer types...")

    for maf_path in maf_files:
        # Пытаемся найти Case ID в пути к файлу или в метаданных по File ID
        file_id = os.path.basename(os.path.dirname(maf_path))
        case_id = meta[meta['File ID'] == file_id]['Case ID'].values[0] if file_id in meta['File ID'].values else None
        
        if not case_id: continue
        
        cancer_type = case_to_cancer.get(case_id, "Unknown")

        try:
            with gzip.open(maf_path, 'rt') as f:
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    subset['Real_Cancer_Type'] = cancer_type
                    results.append(subset)
        except:
            continue

    if not results:
        print("No mutations found.")
        return

    final_df = pd.concat(results)
    
    # Оставляем только те, что мы опознали
    final_df = final_df[final_df['Real_Cancer_Type'] != 'Unknown']

    print("\n" + "="*60)
    print("FINAL ACCURATE MUTATION REPORT")
    print("="*60)
    
    # Сводка по генам и РЕАЛЬНЫМ типам рака
    summary = pd.crosstab(final_df['Real_Cancer_Type'], final_df['Hugo_Symbol'])
    print("\n[Table 1] Mutations per Cancer (Identified by Case ID):")
    print(summary)
    
    print("\n" + "-"*60)
    
    # Только тяжелые мутации
    lof_variants = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = final_df[final_df['Variant_Classification'].isin(lof_variants)]
    
    print(f"\n[Table 2] Loss-of-Function Mutations (The 'Gene Killers'):")
    summary_lof = pd.crosstab(dangerous['Real_Cancer_Type'], dangerous['Hugo_Symbol'])
    print(summary_lof)
    
    print("\n" + "="*60)
    final_df.to_csv('results/hla_mutations_final_v4.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
