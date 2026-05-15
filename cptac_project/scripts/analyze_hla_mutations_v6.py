import pandas as pd
import os
import gzip
import glob

def main():
    maf_files = glob.glob(os.path.join('data/', '**/*.maf.gz'), recursive=True)
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    all_results = []
    print(f"Scanning {len(maf_files)} files... Reading barcodes directly from INSIDE files.")

    for path in maf_files:
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
                    # Вытаскиваем тип рака ПРЯМО из баркода пациента
                    def detect_from_barcode(barcode):
                        b = str(barcode).upper()
                        if 'OV' in b: return 'Ovarian_Cancer'
                        if 'BR' in b: return 'Breast_Cancer'
                        if 'CO' in b: return 'Colon_Cancer'
                        return 'Other'
                    
                    subset['Detected_Cancer'] = subset['Tumor_Sample_Barcode'].apply(detect_from_barcode)
                    all_results.append(subset)
        except:
            continue

    if not all_results:
        print("No mutations found at all.")
        return

    final = pd.concat(all_results)
    
    print("\n" + "="*60)
    print("FINAL SUCCESSFUL REPORT (Direct Barcode Detection)")
    print("="*60)
    
    # Сводка 1: Все мутации
    summary = pd.crosstab(final['Detected_Cancer'], final['Hugo_Symbol'])
    print("\n[Table 1] All HLA/B2M Mutations:")
    print(summary)
    
    print("\n" + "-"*60)
    
    # Сводка 2: Тяжелые поломки (LoF)
    lof = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Splice_Site']
    dangerous = final[final['Variant_Classification'].isin(lof)]
    print("\n[Table 2] SEVERE Loss-of-Function Mutations (The Real Killers):")
    print(pd.crosstab(dangerous['Detected_Cancer'], dangerous['Hugo_Symbol']))
    
    final.to_csv('results/hla_mutations_final_v6.tsv', sep='\t', index=False)
    print("\n" + "="*60)

if __name__ == "__main__":
    main()
