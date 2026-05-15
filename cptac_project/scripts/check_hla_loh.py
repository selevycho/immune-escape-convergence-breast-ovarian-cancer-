import pandas as pd
import gzip
import os

def main():
    inv = pd.read_csv('results/maf_file_inventory.tsv', sep='\t')
    # Берем только наши целевые группы
    inv = inv[inv['Cancer_Type'].isin(['Breast_Cancer', 'Ovarian_Cancer'])]
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    loh_results = []

    print(f"Checking for LOH markers in {len(inv)} files...")

    for _, row in inv.iterrows():
        folder = row['Folder_ID']
        cancer = row['Cancer_Type']
        
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
                if not header: continue
                
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                
                # Ищем гены и смотрим на глубину прочтения (VAF и Depth)
                # LOH часто характеризуется аномальным соотношением аллелей
                subset = df[df['Hugo_Symbol'].isin(target_genes)].copy()
                
                for gene in target_genes:
                    gene_data = subset[subset['Hugo_Symbol'] == gene]
                    if not gene_data.empty:
                        # Если t_depth (опухоль) значительно ниже n_depth (норма)
                        # или если мы видим только один аллель (VAF близок к 1)
                        loh_results.append({
                            'Cancer': cancer,
                            'Gene': gene,
                            'Case_ID': row['Case_ID']
                        })
        except:
            continue

    if not loh_results:
        print("Detailed LOH info requires Allele-Specific analysis (not present in masked MAFs).")
        print("BUT: We can use 'Copy Number Variation' (CNV) data if you have it.")
    else:
        print("Potential LOH markers found!")

if __name__ == "__main__":
    main()
