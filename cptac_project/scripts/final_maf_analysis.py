import pandas as pd
import gzip
import os

def main():
    # 1. Загружаем наш инвентарь
    inv = pd.read_csv('results/maf_file_inventory.tsv', sep='\t')
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    lof_variants = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    
    all_muts = []
    
    print("Processing files based on inventory...")
    
    for _, row in inv.iterrows():
        folder = row['Folder_ID']
        cancer = row['Cancer_Type']
        case = row['Case_ID']
        
        # Ищем файл внутри папки
        data_path = f"data/{folder}/"
        files = [f for f in os.listdir(data_path) if f.endswith('.maf.gz')]
        if not files: continue
        
        maf_path = os.path.join(data_path, files[0])
        
        try:
            with gzip.open(maf_path, 'rt') as f:
                header = None
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                if not header: continue
                
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                # Фильтруем наши гены
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    subset['Cancer_Type'] = cancer
                    subset['Case_ID'] = case
                    # Помечаем, серьезная ли мутация (LoF)
                    subset['Is_LoF'] = subset['Variant_Classification'].isin(lof_variants)
                    all_results.append(subset)
        except:
            continue

    if not all_results:
        print("No mutations found in target genes.")
        return

    final_df = pd.concat(all_results)
    
    print("\n" + "="*50)
    print("FINAL SCIENTIFIC MAF REPORT")
    print("="*50)
    
    # Сводная таблица по типам рака и генам (все мутации)
    summary = pd.crosstab(final_df['Cancer_Type'], final_df['Hugo_Symbol'])
    print("\n[Table 1] Total mutation count:")
    print(summary)
    
    # Сводная таблица только по тяжелым мутациям (LoF)
    lof_df = final_df[final_df['Is_LoF'] == True]
    summary_lof = pd.crosstab(lof_df['Cancer_Type'], lof_df['Hugo_Symbol'])
    print("\n[Table 2] Severe LoF mutations (Gene Killers):")
    print(summary_lof)

    # Расчет частоты (%)
    print("\n[Table 3] Mutation Frequency (% of patients affected):")
    counts = inv['Cancer_Type'].value_counts()
    for cancer in summary.index:
        total_cases = counts[cancer]
        mut_cases = final_df[final_df['Cancer_Type'] == cancer]['Case_ID'].nunique()
        print(f"{cancer}: {round(mut_cases/total_cases*100, 2)}% patients have HLA/B2M mutations")

    final_df.to_csv('results/final_mutation_dataset.tsv', sep='\t', index=False)

if __name__ == "__main__":
    # Исправляем маленькую опечатку в коде (переменная all_results)
    all_results = []
    main()
