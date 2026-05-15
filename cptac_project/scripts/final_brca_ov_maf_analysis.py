import pandas as pd
import gzip
import os

def main():
    # 1. Загружаем инвентарь
    inv = pd.read_csv('results/maf_file_inventory.tsv', sep='\t')
    
    # 2. Оставляем ТОЛЬКО Breast и Ovarian
    target_cancer = ['Breast_Cancer', 'Ovarian_Cancer']
    inv_filtered = inv[inv['Cancer_Type'].isin(target_cancer)].copy()
    
    print(f"Starting analysis for {len(inv_filtered)} targeted files (Breast & Ovarian only)...")
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    lof_variants = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    
    all_results = []

    for _, row in inv_filtered.iterrows():
        folder = row['Folder_ID']
        cancer = row['Cancer_Type']
        case = row['Case_ID']
        
        data_path = f"data/{folder}/"
        try:
            files = [f for f in os.listdir(data_path) if f.endswith('.maf.gz')]
            if not files: continue
            
            maf_path = os.path.join(data_path, files[0])
            
            with gzip.open(maf_path, 'rt') as f:
                header = None
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                if not header: continue
                
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                
                # Ищем наши гены
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    subset['Cancer_Type'] = cancer
                    subset['Case_ID'] = case
                    subset['Is_LoF'] = subset['Variant_Classification'].isin(lof_variants)
                    all_results.append(subset)
        except Exception as e:
            continue

    print("\n" + "="*60)
    print("CLEAN REPORT: BREAST & OVARIAN CANCER ONLY")
    print("="*60)

    if not all_results:
        print("No mutations found in HLA/B2M for these cancer types.")
    else:
        final_df = pd.concat(all_results)
        
        # Таблица 1: Все мутации
        print("\n[Table 1] All detected mutations:")
        summary = pd.crosstab(final_df['Cancer_Type'], final_df['Hugo_Symbol'])
        print(summary)
        
        # Таблица 2: Тяжелые (LoF)
        print("\n[Table 2] Severe LoF Mutations (Gene Killers):")
        lof_df = final_df[final_df['Is_LoF'] == True]
        if not lof_df.empty:
            summary_lof = pd.crosstab(lof_df['Cancer_Type'], lof_df['Hugo_Symbol'])
            print(summary_lof)
        else:
            print("No LoF mutations found.")

        # Таблица 3: Статистика по пациентам
        print("\n[Table 3] Clinical Summary:")
        for cancer in target_cancer:
            total_patients = len(inv_filtered[inv_filtered['Cancer_Type'] == cancer])
            mutated_cases = final_df[final_df['Cancer_Type'] == cancer]['Case_ID'].nunique()
            percent = round((mutated_cases / total_patients) * 100, 2)
            print(f"{cancer}: {mutated_cases} out of {total_patients} patients have mutations ({percent}%)")

        final_df.to_csv('results/final_brca_ov_mutations.tsv', sep='\t', index=False)
        print("\nResults saved to results/final_brca_ov_mutations.tsv")

if __name__ == "__main__":
    main()
