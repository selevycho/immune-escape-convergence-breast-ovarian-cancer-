import pandas as pd
import os

# 1. Setup absolute paths for the cluster environment
HOME = os.path.expanduser("~")
CPTAC_PATH = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")
AURORA_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")

def compare_stealth_genes():
    # Verify file existence
    if not os.path.exists(CPTAC_PATH) or not os.path.exists(AURORA_PATH):
        print("[-] Error: One of the files is missing. Please check the paths.")
        return

    print("[+] Loading CPTAC and AURORA high silencing datasets...")
    
    # 2. Load CPTAC using the specific column names found in diagnostic
    cptac = pd.read_csv(CPTAC_PATH)
    # Both files are already filtered for 'High Silencing' based on previous steps
    # We extract the unique sets of gene names
    cptac_high_genes = set(cptac['Hugo_Symbol'].unique())

    # 3. Load AURORA using its specific column names
    aurora = pd.read_csv(AURORA_PATH)
    aurora_high_genes = set(aurora['gene'].unique())

    # 4. Find the intersection (genes that are hidden in both primary and metastatic)
    shared_hidden = cptac_high_genes.intersection(aurora_high_genes)

    print("\n" + "="*60)
    print("LINE 1 COMPARISON: SHARED HIGH SILENCING GENES")
    print("="*60)
    print(f"Total CPTAC High Silencing genes:   {len(cptac_high_genes)}")
    print(f"Total AURORA High Silencing genes:  {len(aurora_high_genes)}")
    print("-" * 60)
    print(f"SHARED 'STEALTH' GENES FOUND:       {len(shared_hidden)}")
    print("-" * 60)

    # 5. Display a sample of shared genes for verification
    if shared_hidden:
        print("\n[!] Top 20 Shared Hidden Genes:")
        print(sorted(list(shared_hidden))[:20])

    # 6. Save the results for the Jupyter Notebook / Dissertation
    results_dir = os.path.dirname(AURORA_PATH)
    output_file = os.path.join(results_dir, "shared_high_silencing_genes.csv")
    pd.DataFrame({'gene': sorted(list(shared_hidden))}).to_csv(output_file, index=False)
    
    print(f"\n[+] Comparison complete. Shared list saved to: {output_file}")

if __name__ == "__main__":
    compare_stealth_genes()
