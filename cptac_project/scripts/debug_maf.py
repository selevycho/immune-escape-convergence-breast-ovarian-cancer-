import pandas as pd
import gzip
import glob
import os

def main():
    maf_files = glob.glob(os.path.join('data/', '**/*.maf.gz'), recursive=True)
    if not maf_files:
        print("No MAF files found!")
        return

    # Берем первый попавшийся файл для теста
    test_file = maf_files[0]
    print(f"Testing file: {test_file}")

    with gzip.open(test_file, 'rt') as f:
        for line in f:
            if line.startswith('Hugo_Symbol'):
                header = line.strip().split('\t')
                break
        df = pd.read_csv(f, sep='\t', names=header, comment='#', nrows=5)
        
    print("\n--- FIRST 5 ROWS IN MAF ---")
    print(df[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification']])
    
    print("\n--- WHAT WE HAVE IN METADATA ---")
    meta = pd.read_csv('manifests/metadata_final.tsv', sep='\t', nrows=5)
    print(meta[['Case ID', 'Project ID']])

if __name__ == "__main__":
    main()
