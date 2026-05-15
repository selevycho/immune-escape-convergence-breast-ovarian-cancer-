import pandas as pd
import os
import sys

def main():
    # 1. Path setup
    metadata_path = 'manifests/metadata_final.tsv'
    data_dir = 'data/'
    output_dir = 'results'
    output_path = os.path.join(output_dir, 'hla_expression_matrix.tsv')

    print("--- Starting HLA Expression Extraction ---")

    # 2. Safety checks
    if not os.path.exists(metadata_path):
        print(f"Error: Metadata file {metadata_path} not found!")
        sys.exit(1)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 3. Load metadata and filter for RNA-Seq files only
    metadata = pd.read_csv(metadata_path, sep='\t')
    
    # We only need 'Gene Expression Quantification' files (the .tsv ones)
    rna_files = metadata[metadata['Data Type'] == 'Gene Expression Quantification'].copy()
    
    if rna_files.empty:
        print("Error: No RNA-Seq (Gene Expression) files found in metadata!")
        sys.exit(1)

    print(f"Found {len(rna_files)} expression files. Starting extraction...")

    extracted_data = []
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']

    # 4. Processing each file
    for idx, row in rna_files.iterrows():
        file_id = row['File ID']
        file_name = row['File Name']
        
        # Path to the actual data file
        # Note: GDC stores files as data/<file_id>/<file_name>
        full_path = os.path.join(data_dir, file_id, file_name)
        
        if not os.path.exists(full_path):
            # Try to see if it's just in data/<file_id>/
            # Sometimes gdc-client behavior varies
            continue

        try:
            # GDC TSV files have 1-5 lines of headers/metadata
            # We look for the line where gene_id or gene_name starts
            # Usually skipping 1 or 6 lines works
            df_gene = pd.read_csv(full_path, sep='\t', skiprows=1)
            
            # If the column name is wrong due to GDC version, we adjust
            if 'gene_name' not in df_gene.columns:
                df_gene = pd.read_csv(full_path, sep='\t', skiprows=6)

            # Extract TPM (normalized expression) for our genes
            # Column 'tpm_unstranded' is best for comparing between samples
            subset = df_gene[df_gene['gene_name'].isin(target_genes)]
            
            # Create a dictionary for this sample with basic info
            sample_info = {
                'Case_ID': row['Case ID'],
                'Project': row['Project ID'],
                'Tissue_Type': row['Tissue Type'],
                'File_ID': file_id
            }
            
            # Add expression values to our dictionary
            for _, g_row in subset.iterrows():
                gene_name = g_row['gene_name']
                tpm_value = g_row['tpm_unstranded']
                sample_info[gene_name] = tpm_value
            
            # Check if we actually found all genes
            if all(gene in sample_info for gene in target_genes):
                extracted_data.append(sample_info)
            
        except Exception as e:
            print(f"Warning: Could not process file {file_id}: {e}")

    # 5. Finalize and Save
    if not extracted_data:
        print("Error: No gene data could be extracted. Check file formats!")
        sys.exit(1)

    final_df = pd.DataFrame(extracted_data)
    
    # Save results
    final_df.to_csv(output_path, sep='\t', index=False)
    
    print("\n" + "="*50)
    print(f"SUCCESS: Extracted data for {len(final_df)} samples.")
    print(f"Matrix saved to: {output_path}")
    print("="*50)
    print("\nPreview of extracted HLA levels:")
    print(final_df.head())

if __name__ == "__main__":
    main()
