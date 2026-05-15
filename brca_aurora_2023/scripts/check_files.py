import pandas as pd
import os

# Define absolute paths
HOME = os.path.expanduser("~")
CPTAC_FILE = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")
AURORA_FILE = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")

def check_file(path, name):
    print(f"\n{'='*20} CHECKING: {name} {'='*20}")
    if not os.path.exists(path):
        print(f"❌ Error: File not found at {path}")
        return

    # Load file with automatic delimiter detection
    try:
        df = pd.read_csv(path, sep=None, engine='python')
        print(f"📍 Path: {path}")
        print(f"📊 Shape: {df.shape} (rows, columns)")
        print(f"📝 Column names: {list(df.columns)}")
        print("\n--- TOP 10 ROWS ---")
        print(df.head(10).to_string(index=False))
    except Exception as e:
        print(f"❌ Failed to read file: {e}")

if __name__ == "__main__":
    # Check Project 1 (Primary)
    check_file(CPTAC_FILE, "CPTAC (Primary)")
    
    # Check Project 2 (Metastatic)
    check_file(AURORA_FILE, "AURORA (Metastatic)")
