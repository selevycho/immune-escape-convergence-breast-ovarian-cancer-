import pandas as pd
import os

# 1. Setup absolute paths
HOME = os.path.expanduser("~")
# CPTAC (Primary) - Using the path and columns from your screenshot
CPTAC_FILE = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")
# AURORA (Metastatic) - Results from our current work
AURORA_FILE = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")

def compare_high_silencing():
    # Verify both files exist
    if not os.path.exists(CPTAC_FILE):
        print(f"[-] Error: CPTAC file not found at {CPTAC_FILE}")
        return
    if not os.path.exists(AURORA_FILE):
        print(f"[-] Error: Aurora file not found at {AURORA_FILE}")
        return

    print("[+] Loading and mapping datasets...")
    
    # 2. Load CPTAC and map columns:
    # hugo_symbol -> gene
    # silenced_percent -> ratio
    cptac = pd.read_csv(CPTAC_FILE)
    
    # Filter for high silencing (> 50%) using your specific column name 'silenced_percent'
    # and identifying genes via 'hugo_symbol'
    cptac_high_genes = set(cptac[cptac['silenced_percent'] > 50]['hugo_symbol'])

    # 3. Load AURORA
    # Here the columns are already named 'gene' and it's already filtered for high silencing
    aurora = pd.read_csv(AURORA_FILE)
    aurora_high_genes = set(aurora['gene'])

    # 4. Perform Intersection
    shared_hidden = cptac_high_genes.intersection(aurora_high_genes)

    print("\n" + "="*60)
    print("FINAL COMPARISON: HIGH SILENCING (STEALTH MODE)")
    print("="*60)
    print(f"CPTAC (Primary) High Silencing genes:   {len(cptac_high_genes)}")
    print(f"AURORA (Metastatic) High Silencing genes: {len(aurora_high_genes)}")
    print("-" * 60)
    print(f"SHARED HIGH SILENCING GENES:           {len(shared_hidden)}")
    print("-" * 60)

    # Display results
    if shared_hidden:
        print("\n[!] Top shared 'Stealth' genes found in both datasets (first 30):")
        print(sorted(list(shared_hidden))[:30])

    # 5. Save the intersection list for your Jupyter Notebook
    output_path = os.path.join(os.path.dirname(AURORA_FILE), "shared_high_silencing_genes.csv")
    pd.DataFrame({'gene': sorted(list(shared_hidden))}).to_csv(output_path, index=False)
    print(f"\n[+] Results successfully saved to: {output_path}")

if __name__ == "__main__":
    compare_high_silencing()
