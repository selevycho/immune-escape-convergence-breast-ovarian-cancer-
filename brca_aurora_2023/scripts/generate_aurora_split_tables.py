import pandas as pd
import os

# 1. Setup absolute paths based on the script's location
# This ensures it works regardless of where you launch it from
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "results"))
INPUT_FILE = os.path.join(RESULTS_DIR, "line1_aurora_merged_data.csv")

def generate_split_tables():
    # Check if the required merged file exists from the previous step
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Error: {INPUT_FILE} not found. Please ensure you ran the Line 1 merging script first.")
        return

    print("🚀 Loading processed AURORA mutation data...")
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Define silencing threshold (Z < -1.0 is our standard for 'Low Expression')
    # A mutation is considered 'silenced' if its gene expression is significantly below average
    threshold = -1.0
    df['is_silenced'] = df['Z_score'] < threshold

    # 3. Aggregate data per gene to calculate silencing rates
    # We need: total mutations per gene, count of silenced mutations, and the ratio (%)
    full_table = df.groupby('Hugo_Symbol').agg(
        N_mutations=('is_silenced', 'count'),
        N_silenced=('is_silenced', 'sum')
    ).reset_index()
    
    full_table.rename(columns={'Hugo_Symbol': 'gene'}, inplace=True)
    full_table['ratio'] = (full_table['N_silenced'] / full_table['N_mutations']) * 100

    # 4. Split into High and Low Silencing tables
    # High Silencing: Genes where more than 50% of mutations are hidden (Ratio > 50%)
    high_silencing_df = full_table[full_table['ratio'] > 50].copy()
    
    # Low Silencing (Exposed): Genes where mutations are mostly visible (Ratio < 20%)
    # This often includes critical driver genes like TP53 or PIK3CA
    low_silencing_df = full_table[full_table['ratio'] < 20].copy()

    # Sort for better readability (by mutation count then ratio)
    high_silencing_df = high_silencing_df.sort_values(by=['N_mutations', 'ratio'], ascending=False)
    low_silencing_df = low_silencing_df.sort_values(by=['N_mutations', 'ratio'], ascending=True)

    # 5. Save the final CSV files to the results folder
    high_path = os.path.join(RESULTS_DIR, "high_silencing_aurora_all.csv")
    low_path = os.path.join(RESULTS_DIR, "low_silencing_aurora_all.csv")

    high_silencing_df.to_csv(high_path, index=False)
    low_silencing_df.to_csv(low_path, index=False)
    
    print("-" * 50)
    print(f"✅ Success! Two tables created in: {RESULTS_DIR}")
    print(f"1. High Silencing (Hidden): {len(high_silencing_df)} genes")
    print(f"2. Low Silencing (Exposed): {len(low_silencing_df)} genes")
    print("-" * 50)

if __name__ == "__main__":
    generate_split_tables()
