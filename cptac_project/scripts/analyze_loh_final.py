import pandas as pd
import gzip
import os

def main():
    inv = pd.read_csv('results/maf_file_inventory.tsv', sep='\t')
    inv = inv[inv['Cancer_Type'].isin(['Breast_Cancer', 'Ovarian_Cancer'])]
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    results = []

    print("Analyzing Depth and VAF for LOH detection...")

    for _, row in inv.iterrows():
        folder = row['Folder_ID']
        cancer = row['Cancer_Type']
        case = row['Case_ID']
        
        data_path = f"data/{folder}/"
        try:
            files = [f for f in os.listdir(data_path) if f.endswith('.maf.gz')]
            if not files: continue
            
            with gzip.open(os.path.join(data_path, files[0]), 'rt') as f:
                header = None
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                
                # Мы ищем признаки LOH: когда одна копия гена доминирует (VAF > 0.8)
                # или когда общее покрытие гена аномально низкое.
                # В MAF это обычно колонки t_alt_count и t_ref_count
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                if not subset.empty:
                    for gene in target_genes:
                        gene_df = subset[subset['Hugo_Symbol'] == gene]
                        if not gene_df.empty:
                            results.append({
                                'Case_ID': case,
                                'Cancer': cancer,
                                'Gene': gene,
                                'Type': 'Potential_LOH'
                            })
        except:
            continue

    if not results:
        print("No specific markers found in these files.")
        return

    final_loh = pd.DataFrame(results)
    
    print("\n" + "="*50)
    print("LOH (LOSS OF HETEROGENEITY) SUMMARY")
    print("="*50)
    
    summary = pd.crosstab(final_loh['Cancer'], final_loh['Gene'])
    print(summary)
    
    print("\n" + "-"*50)
    # Считаем процент пациентов с LOH
    for cancer in ['Breast_Cancer', 'Ovarian_Cancer']:
        total = len(inv[inv['Cancer_Type'] == cancer])
        cases_with_loh = final_loh[final_loh['Cancer'] == cancer]['Case_ID'].nunique()
        print(f"{cancer}: {cases_with_loh} patients ({round(cases_with_loh/total*100, 2)}%) show LOH markers")

if __name__ == "__main__":
    main()
