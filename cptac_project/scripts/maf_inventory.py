import pandas as pd
import os
import glob

def main():
    meta_path = 'manifests/metadata_final.tsv'
    if not os.path.exists(meta_path):
        print(f"Error: {meta_path} not found!")
        return
        
    meta = pd.read_csv(meta_path, sep='\t')
    
    # ИСПРАВЛЕНО: добавляем .str перед .strip()
    file_to_project = dict(zip(meta['File ID'].str.strip(), meta['Project ID'].str.strip()))
    file_to_case = dict(zip(meta['File ID'].str.strip(), meta['Case ID'].str.strip()))

    maf_paths = glob.glob(os.path.join('data/', '**/*.maf.gz'), recursive=True)
    
    inventory = []
    print(f"Starting inventory of {len(maf_paths)} files...\n")

    for path in maf_paths:
        folder_id = os.path.basename(os.path.dirname(path)).strip()
        
        project = file_to_project.get(folder_id, "NOT_FOUND_IN_METADATA")
        case_id = file_to_case.get(folder_id, "UNKNOWN")
        
        cancer_type = "Other"
        if 'BR' in str(case_id): cancer_type = 'Breast_Cancer'
        elif 'OV' in str(case_id): cancer_type = 'Ovarian_Cancer'
        elif 'CO' in str(case_id): cancer_type = 'Colon_Cancer'
        
        inventory.append({
            'Folder_ID': folder_id,
            'Project_ID': project,
            'Case_ID': case_id,
            'Cancer_Type': cancer_type
        })

    df_inv = pd.DataFrame(inventory)
    
    print("--- INVENTORY SUMMARY ---")
    summary = df_inv.groupby(['Project_ID', 'Cancer_Type']).size().reset_index(name='File_Count')
    print(summary)
    
    print(f"\nFiles not found in metadata: {len(df_inv[df_inv['Project_ID'] == 'NOT_FOUND_IN_METADATA'])}")
    
    df_inv.to_csv('results/maf_file_inventory.tsv', sep='\t', index=False)
    print("\nInventory saved to results/maf_file_inventory.tsv")

if __name__ == "__main__":
    main()
