import pandas as pd
import os
import sys

def main():
    sample_sheet_path = 'manifests/gdc_sample_sheet.2026-05-07.tsv'
    data_dir = 'data/'
    output_path = 'manifests/metadata_final.tsv'

    print("--- Metadata Processing Started ---")

    # Load Sample Sheet
    if not os.path.exists(sample_sheet_path):
        print(f"Error: {sample_sheet_path} not found!")
        sys.exit(1)
    
    df = pd.read_csv(sample_sheet_path, sep='\t')

    # Identify downloaded folders
    downloaded_ids = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d))]
    
    # Filter the dataframe
    df_final = df[df['File ID'].isin(downloaded_ids)].copy()

    if df_final.empty:
        print("Error: No matches found. Check your data directory!")
        sys.exit(1)

    # Use the specific columns we found in your output
    # Project ID and Tissue Type are the keys here
    print("\n" + "="*50)
    print("DOWNLOADED DATA STATISTICS:")
    print("="*50)
    
    # We group by Project and Tissue Type (Tumor vs Normal)
    summary = df_final.groupby(['Project ID', 'Tissue Type']).size()
    print(summary)
    
    print("="*50)

    # Save the cleaned metadata
    df_final.to_csv(output_path, sep='\t', index=False)
    print(f"SUCCESS: Final metadata saved to: {output_path}")

if __name__ == "__main__":
    main()
