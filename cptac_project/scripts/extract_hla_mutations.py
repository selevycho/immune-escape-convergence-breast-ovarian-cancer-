import pandas as pd
import os
import gzip
import glob

def main():
    data_dir = 'data/'
    output_path = 'results/hla_mutations_raw.tsv'
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    # Ищем все файлы .maf.gz в подпапках
    maf_files = glob.glob(os.path.join(data_dir, '**/*.maf.gz'), recursive=True)
    
    print(f"Found {len(maf_files)} MAF files. Starting scan...")
    
    all_mutations = []

    for maf_path in maf_files:
        try:
            # Читаем сжатый файл, пропускаем комментарии (строки на #)
            with gzip.open(maf_path, 'rt') as f:
                # Находим строку с заголовками
                for line in f:
                    if line.startswith('Hugo_Symbol'):
                        header = line.strip().split('\t')
                        break
                
                # Читаем данные
                df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                
                # Фильтруем только наши гены
                subset = df[df['Hugo_Symbol'].isin(target_genes)]
                
                if not subset.empty:
                    # Оставляем только важные колонки для анализа
                    cols_to_keep = ['Hugo_Symbol', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode']
                    all_mutations.append(subset[cols_to_keep])
                    
        except Exception as e:
            print(f"Error processing {maf_path}: {e}")

    if all_mutations:
        final_df = pd.concat(all_mutations)
        final_df.to_csv(output_path, sep='\t', index=False)
        print(f"\nSUCCESS: Found {len(final_df)} mutations in target genes.")
        print(final_df['Hugo_Symbol'].value_counts())
    else:
        print("No mutations found in HLA genes or B2M.")

if __name__ == "__main__":
    main()
