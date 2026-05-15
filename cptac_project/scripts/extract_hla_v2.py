import pandas as pd
import os

def main():
    metadata = pd.read_csv('manifests/metadata_final.tsv', sep='\t')
    data_dir = 'data/'
    
    # Filter only RNA-seq files
    rna = metadata[metadata['Data Type'] == 'Gene Expression Quantification'].copy()
    
    results = []
    genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    print(f"Processing {len(rna)} RNA-seq files...")
    for _, row in rna.iterrows():
        path = os.path.join(data_dir, row['File ID'], row['File Name'])
        if os.path.exists(path):
            try:
                # GDC format check
                df = pd.read_csv(path, sep='\t', skiprows=1)
                if 'gene_name' not in df.columns:
                    df = pd.read_csv(path, sep='\t', skiprows=6)
                
                subset = df[df['gene_name'].isin(genes)]
                res = {
                    'Case_ID': row['Case ID'],
                    'Project': row['Project ID'],
                    'Tissue': row['Tissue Type']
                }
                for _, g in subset.iterrows():
                    res[g['gene_name']] = g['tpm_unstranded']
                results.append(res)
            except:
                continue

    final_df = pd.DataFrame(results)
    final_df.to_csv('results/hla_matrix_v2.tsv', sep='\t', index=False)
    print("Success: Matrix v2 saved to results/hla_matrix_v2.tsv")

if __name__ == "__main__":
    main()
