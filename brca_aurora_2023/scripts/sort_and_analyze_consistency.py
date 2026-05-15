import pandas as pd
import os

# 1. Setup paths
HOME = os.path.expanduser("~")
INPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/combined_low_silencing_stats.csv")
OUTPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/sorted_stability_exposed_genes.csv")

def refine_and_sort():
    if not os.path.exists(INPUT_PATH):
        print(f"[-] Error: Input file not found at {INPUT_PATH}")
        return

    # Load the combined data
    df = pd.read_csv(INPUT_PATH)

    # 2. Calculate New Metrics for Visualization
    # 'Ratio_Diff' shows how much the silencing behavior changed during metastasis
    df['Ratio_Diff'] = abs(df['CPTAC_Percent'] - df['AURORA_Percent'])
    
    # 'Mean_Exposure' shows the average silencing percentage (lower is more "naked")
    df['Mean_Percent'] = (df['CPTAC_Percent'] + df['AURORA_Percent']) / 2
    
    # 'Total_Impact' combines mutation counts from both cohorts
    df['Total_Impact'] = df['CPTAC_Mutations'] + df['AURORA_Mutations']

    # 3. Categorize genes based on their consistency
    # We define 'Highly Consistent' as genes with less than 5% difference in silencing
    def categorize(row):
        if row['Ratio_Diff'] < 5 and row['Mean_Percent'] < 5:
            return "1. Highly Consistent & Exposed"
        elif row['Ratio_Diff'] < 10:
            return "2. Stable Target"
        else:
            return "3. Variable Exposure"

    df['Consistency_Category'] = df.apply(categorize, axis=1)

    # 4. Final Sorting
    # We want: 
    # - Consistency category first
    # - Lowest Mean_Percent second (the most "naked" genes)
    # - Highest Total_Impact third (the most frequent mutations)
    df_sorted = df.sort_values(
        by=['Consistency_Category', 'Mean_Percent', 'Total_Impact'],
        ascending=[True, True, False]
    )

    # 5. Reorder columns for better readability
    cols = [
        'gene', 'Consistency_Category', 'Mean_Percent', 'Ratio_Diff', 
        'CPTAC_Percent', 'AURORA_Percent', 'CPTAC_Mutations', 'AURORA_Mutations'
    ]
    df_sorted = df_sorted[cols]

    # 6. Save and Print
    df_sorted.to_csv(OUTPUT_PATH, index=False)
    
    print("\n" + "="*80)
    print("REFINED GENE STABILITY ANALYSIS (LINE 1)")
    print("="*80)
    print(f"{'Gene':<12} | {'Category':<30} | {'Mean %':<8} | {'Diff':<6}")
    print("-" * 80)
    print(df_sorted.head(20).to_string(index=False, header=False))
    print("-" * 80)
    print(f"[+] Beautifully sorted table saved to: {OUTPUT_PATH}")

if __name__ == "__main__":
    refine_and_sort()
