import pandas as pd
import os

# 1. Setup absolute paths based on the script location
# This prevents PermissionError when running from different directories
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "results"))

# Create results directory if it doesn't exist
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

def run_analysis():
    print("🚀 Loading AURORA 2023 datasets...")
    
    mut_file = os.path.join(DATA_DIR, "data_mutations.txt")
    rna_file = os.path.join(DATA_DIR, "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")

    # Check if files exist before processing
    if not os.path.exists(mut_file) or not os.path.exists(rna_file):
        print(f"❌ Error: Files not found in {DATA_DIR}")
        return

    # Load mutation data
    # We only need Hugo_Symbol and Tumor_Sample_Barcode
    muts = pd.read_csv(mut_file, sep='\t', low_memory=False, 
                       usecols=['Hugo_Symbol', 'Tumor_Sample_Barcode'])
    
    # Load RNA Z-scores
    # Rows = Genes, Columns = Samples
    rna = pd.read_csv(rna_file, sep='\t')
    
    print("🧬 Processing and merging data...")
    
    # Melt RNA dataframe to convert it to long format for easier merging
    # Resulting columns: Hugo_Symbol, Sample_ID, Z_score
    rna_long = rna.melt(id_vars=['Hugo_Symbol'], var_name='Sample_ID', value_name='Z_score')
    
    # Perform inner merge to match mutations with their respective expression levels
    # We merge on Gene name and Sample ID
    merged = pd.merge(
        muts, 
        rna_long, 
        left_on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], 
        right_on=['Hugo_Symbol', 'Sample_ID'], 
        how='inner'
    )

    # Drop missing values (mutations without corresponding RNA data)
    merged = merged.dropna(subset=['Z_score'])
    total_valid = len(merged)

    print(f"\n--- LINE 1 RESULTS: AURORA 2023 (N={total_valid} mutations) ---")
    
    # 2. Calculate Silencing Ratios based on Z-score thresholds
    # Z < -1.0: Low expression (Bottom 16% of distribution)
    # Z < -2.0: Deep silencing (Bottom 2.3% of distribution)
    thresholds = {
        'Low Expression (Z < -1.0)': -1.0,
        'Deep Silencing (Z < -2.0)': -2.0
    }

    summary_stats = []
    for label, t in thresholds.items():
        silenced_count = len(merged[merged['Z_score'] < t])
        ratio = (silenced_count / total_valid) * 100
        print(f"{label}: {silenced_count} muts ({ratio:.2f}%)")
        summary_stats.append({'Threshold': label, 'Value': t, 'Count': silenced_count, 'Ratio': ratio})

    # 3. Save results for further use in Jupyter Notebook
    # Detailed table for plotting
    output_merged = os.path.join(RESULTS_DIR, "line1_aurora_merged_data.csv")
    merged.to_csv(output_merged, index=False)
    
    # Summary table for the final report
    output_summary = os.path.join(RESULTS_DIR, "line1_aurora_summary.csv")
    pd.DataFrame(summary_stats).to_csv(output_summary, index=False)

    print(f"\n✅ Success! Detailed results saved to: {RESULTS_DIR}")

if __name__ == "__main__":
    run_analysis()
