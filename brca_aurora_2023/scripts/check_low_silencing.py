import pandas as pd
import os

# Define absolute paths
HOME = os.path.expanduser("~")

# We expect the file to be named low_silencing_all.csv based on your correction
CPTAC_LOW_PATH = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/low_silencing_all.csv")
AURORA_LOW_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/low_silencing_aurora_all.csv")

def diagnostic_check():
    print(f"\n{'='*30} DIAGNOSTICS: LOW SILENCING (RETRY) {'='*30}")
    
    # 1. Check CPTAC Low Silencing
    if os.path.exists(CPTAC_LOW_PATH):
        print(f"\n[+] File Found: CPTAC Low Silencing")
        df_cptac = pd.read_csv(CPTAC_LOW_PATH, sep=None, engine='python')
        print(f"📍 Columns: {list(df_cptac.columns)}")
        print("--- TOP 10 ROWS ---")
        print(df_cptac.head(10).to_string(index=False))
    else:
        print(f"❌ Error: File NOT FOUND at {CPTAC_LOW_PATH}")
        # Let's list files in that directory to be sure
        cptac_res_dir = os.path.dirname(CPTAC_LOW_PATH)
        if os.path.exists(cptac_res_dir):
            print(f"\n📂 Files actually present in {cptac_res_dir}:")
            print(os.listdir(cptac_res_dir))

    # 2. Check AURORA Low Silencing (should be fine as before)
    if os.path.exists(AURORA_LOW_PATH):
        print(f"\n[+] File Found: AURORA Metastatic Low Silencing")
        df_aurora = pd.read_csv(AURORA_LOW_PATH, sep=None, engine='python')
        print(f"📍 Columns: {list(df_aurora.columns)}")
    else:
        print(f"❌ Error: AURORA file not found.")

if __name__ == "__main__":
    diagnostic_check()
