import pandas as pd
import os

def main():
    metadata = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    data_dir = 'data/'
    
    # Filter only RNA-seq
    rna = metadata[metadata['Data Type'] == 'Gene Expression Quantification'].copy()
    
    # DEBUG: Посмотрим, что реально написано в Case ID
    # В CPTAC-2 рак яичников обычно начинается на '17OV', '04OV' и т.д.
    def identify_cancer(row):
        case = str(row['Case ID'])
        if 'OV' in case: return 'Ovarian_Cancer'
        if 'BR' in case: return 'Breast_Cancer'
        if 'CO' in case: return 'Colon_Cancer'
        return 'Other'

    rna['Real_Project'] = rna.apply(identify_cancer, axis=1)
    
    results = []
    genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    print("Detected types in your data:")
    print(rna['Real_Project'].value_counts())

    for _, row in rna.iterrows():
        path = os.path.join(data_dir, row['File ID'], row['File Name'])
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep='\t', skiprows=1)
                if 'gene_name' not in df.columns:
                    df = pd.read_csv(path, sep='\t', skiprows=6)
                
                subset = df[df['gene_name'].isin(genes)]
                res = {
                    'Case_ID': row['Case ID'],
                    'Project': row['Real_Project'] # Используем нашу новую метку
                }
                for _, g in subset.iterrows():
                    res[g['gene_name']] = g['tpm_unstranded']
                results.append(res)
            except:
                continue

    pd.DataFrame(results).to_csv('results/hla_matrix_v3.tsv', sep='\t', index=False)
    print("Success: Matrix v3 saved.")

if __name__ == "__main__":
    main()
