import pandas as pd
import os

# 1. Setup absolute paths
HOME = os.path.expanduser("~")
CPTAC_PATH = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/low_silencing_all.csv")
AURORA_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/low_silencing_aurora_all.csv")
OUTPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/combined_low_silencing_stats.csv")

def generate_combined_table():
    if not os.path.exists(CPTAC_PATH) or not os.path.exists(AURORA_PATH):
        print("[-] Error: Missing input files.")
        return

    print("[+] Loading and merging datasets side-by-side...")

    # 2. Load CPTAC (Primary)
    cptac = pd.read_csv(CPTAC_PATH)
    # Keeping: Hugo_Symbol, Total_Mutations, Silenced_Count, Expressed_Count, Silenced_Percent
    cptac = cptac[['Hugo_Symbol', 'Total_Mutations', 'Silenced_Count', 'Expressed_Count', 'Silenced_Percent']]
    
    # Add prefix to CPTAC columns to distinguish them
    cptac.columns = ['gene', 'CPTAC_Mutations', 'CPTAC_Silenced', 'CPTAC_Expressed', 'CPTAC_Percent']

    # 3. Load AURORA (Metastatic)
    aurora = pd.read_csv(AURORA_PATH)
    
    # Calculate 'Expressed_Count' for Aurora to match CPTAC format
    aurora['N_expressed'] = aurora['N_mutations'] - aurora['N_silenced']
    
    # Selecting and ordering AURORA columns
    aurora = aurora[['gene', 'N_mutations', 'N_silenced', 'N_expressed', 'ratio']]
    
    # Add prefix to AURORA columns
    aurora.columns = ['gene', 'AURORA_Mutations', 'AURORA_Silenced', 'AURORA_Expressed', 'AURORA_Percent']

    # 4. Merge tables (Inner Join on 'gene')
    # This keeps only genes present in BOTH datasets
    combined_df = pd.merge(cptac, aurora, on='gene', how='inner')

    # 5. Sort by importance (Total mutations in both cohorts)
    combined_df['Total_Combined_Mutations'] = combined_df['CPTAC_Mutations'] + combined_df['AURORA_Mutations']
    combined_df = combined_df.sort_values(by='Total_Combined_Mutations', ascending=False)
    
    # Drop the temporary sorting column
    combined_df = combined_df.drop(columns=['Total_Combined_Mutations'])

    # 6. Save results
    combined_df.to_csv(OUTPUT_PATH, index=False)
    
    print("\n" + "="*60)
    print("SIDE-BY-SIDE LOW SILENCING COMPARISON")
    print("="*60)
    print(f"Total Shared Exposed Genes: {len(combined_df)}")
    print("-" * 60)
    print("TOP 10 SHARED GENES (by mutation density):")
    # Displaying top 10 with formatted view
    print(combined_df.head(10).to_string(index=False))
    print("-" * 60)
    print(f"[+] Final combined table saved to: {OUTPUT_PATH}")

if __name__ == "__main__":
    generate_combined_table()
