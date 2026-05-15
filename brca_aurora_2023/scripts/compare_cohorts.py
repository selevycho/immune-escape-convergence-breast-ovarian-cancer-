import pandas as pd
import os

# Define absolute paths for both projects
# We use os.path.expanduser to ensure the home directory is mapped correctly
HOME = os.path.expanduser("~")

# Path to the PREVIOUS project results (Primary Breast Cancer)
CPTAC_FILE = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")

# Path to the CURRENT project results (Metastatic Breast Cancer)
AURORA_DIR = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results")
high_aurora_file = os.path.join(AURORA_DIR, "high_silencing_aurora_all.csv")
low_aurora_file = os.path.join(AURORA_DIR, "low_silencing_aurora_all.csv")

def run_comparison():
    # Verify that all files exist before processing
    if not os.path.exists(CPTAC_FILE):
        print(f"[-] Error: CPTAC file not found at {CPTAC_FILE}")
        return
    if not os.path.exists(high_aurora_file) or not os.path.exists(low_aurora_file):
        print("[-] Error: Aurora results not found. Run Line 1 scripts first.")
        return
    
    # 1. Load the data
    print("[+] Loading CPTAC and AURORA gene lists...")
    cptac = pd.read_csv(CPTAC_FILE)
    aur_high_df = pd.read_csv(high_aurora_file)
    aur_low_df = pd.read_csv(low_aurora_file)

    # 2. Extract gene sets based on silencing status
    # In CPTAC we filter by the 'silencing_status' column
    cptac_high = set(cptac[cptac['silencing_status'] == 'High']['gene'])
    cptac_low = set(cptac[cptac['silencing_status'] == 'Exposed']['gene'])

    # In AURORA we use the pre-split tables
    aur_high = set(aur_high_df['gene'])
    aur_low = set(aur_low_df['gene'])

    print("\n" + "="*60)
    print("ANALYSIS: PRIMARY (CPTAC) VS METASTATIC (AURORA)")
    print("="*60)

    # 3. Intersection: Shared High Silencing (Universal Stealth)
    shared_high = cptac_high.intersection(aur_high)
    print(f"\n[!] SHARED HIGH SILENCING (Always Hidden):")
    print(f"Common genes: {len(shared_high)}")
    if shared_high:
        print(f"Examples: {list(shared_high)[:15]}")

    # 4. Intersection: Shared Low Silencing (Universal Targets)
    shared_low = cptac_low.intersection(aur_low)
    print(f"\n[!] SHARED LOW SILENCING (Always Exposed):")
    print(f"Common genes: {len(shared_low)}")
    if shared_low:
        print(f"Examples: {list(shared_low)[:15]}")

    # 5. Biological Switch: Hidden in Primary -> Exposed in Metastatic
    # This shows genes that become 'vulnerable' during metastatic expansion
    switchers = cptac_high.intersection(aur_low)
    print(f"\n[!] GENES THAT LOST PROTECTION (Hidden -> Exposed):")
    print(f"Count: {len(switchers)}")
    if switchers:
        print(f"Examples: {list(switchers)[:15]}")

    # 6. Save shared exposed genes for Line 3 Convergence analysis
    output_path = os.path.join(AURORA_DIR, "shared_exposed_genes.csv")
    pd.DataFrame({'gene': list(shared_low)}).to_csv(output_path, index=False)
    
    print(f"\n[+] Results saved to: {output_path}")

if __name__ == "__main__":
    run_comparison()
