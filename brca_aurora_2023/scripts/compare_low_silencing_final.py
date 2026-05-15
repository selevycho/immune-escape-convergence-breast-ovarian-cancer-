import pandas as pd
import os

# 1. Setup absolute paths
HOME = os.path.expanduser("~")
CPTAC_LOW_FILE = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/low_silencing_all.csv")
AURORA_LOW_FILE = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/low_silencing_aurora_all.csv")

def compare_exposed_genes():
    # Verify file existence
    if not os.path.exists(CPTAC_LOW_FILE) or not os.path.exists(AURORA_LOW_FILE):
        print("[-] Error: One of the files is missing. Please check the paths.")
        return

    print("[+] Loading CPTAC and AURORA exposed gene datasets...")
    
    # 2. Load CPTAC (Primary)
    # Based on your diagnostic, columns are: 'Hugo_Symbol', 'Silenced_Percent', etc.
    cptac = pd.read_csv(CPTAC_LOW_FILE)
    cptac_exposed_genes = set(cptac['Hugo_Symbol'].unique())

    # 3. Load AURORA (Metastatic)
    # Based on your diagnostic, columns are: 'gene', 'ratio', etc.
    aurora = pd.read_csv(AURORA_LOW_FILE)
    aurora_exposed_genes = set(aurora['gene'].unique())

    # 4. Find the intersection (Universal Targets)
    # These are genes that are exposed (not silenced) in both primary and metastatic stages
    shared_exposed = cptac_exposed_genes.intersection(aurora_exposed_genes)

    print("\n" + "="*60)
    print("LINE 1 COMPARISON: SHARED EXPOSED GENES (Universal Targets)")
    print("="*60)
    print(f"Total CPTAC Exposed genes:   {len(cptac_exposed_genes)}")
    print(f"Total AURORA Exposed genes:  {len(aurora_exposed_genes)}")
    print("-" * 60)
    print(f"SHARED EXPOSED GENES FOUND:  {len(shared_exposed)}")
    print("-" * 60)

    # 5. Display a sample of shared genes
    if shared_exposed:
        print("\n[!] Top 30 Shared Exposed Genes (Alphabetical):")
        sample = sorted(list(shared_exposed))[:30]
        print(sample)

    # 6. Save the shared list for the final Line 3 (B2M) Convergence argument
    results_dir = os.path.dirname(AURORA_LOW_FILE)
    output_file = os.path.join(results_dir, "shared_exposed_genes_final.csv")
    pd.DataFrame({'gene': sorted(list(shared_exposed))}).to_csv(output_file, index=False)
    
    print(f"\n[+] Success! Shared exposed list saved to: {output_file}")

if __name__ == "__main__":
    compare_exposed_genes()
