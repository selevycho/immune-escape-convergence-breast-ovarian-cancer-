import pandas as pd
import os

# 1. Paths
HOME = os.path.expanduser("~")
INPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/combined_low_silencing_stats.csv")
OUTPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/final_line1_stability_report.csv")

def generate_clean_report():
    if not os.path.exists(INPUT_PATH):
        print(f"[-] Error: {INPUT_PATH} not found.")
        return

    df = pd.read_csv(INPUT_PATH)

    # 2. Add metrics
    df['Mean_Pct'] = (df['CPTAC_Percent'] + df['AURORA_Percent']) / 2
    df['Diff'] = abs(df['CPTAC_Percent'] - df['AURORA_Percent'])
    df['Total_Muts'] = df['CPTAC_Mutations'] + df['AURORA_Mutations']

    # 3. Sorting: Prioritize absolute 0% silencing in both cohorts
    # These are the most "dangerous" genes for the tumor
    df_sorted = df.sort_values(by=['Mean_Pct', 'Diff', 'Total_Muts'], ascending=[True, True, False])

    # 4. Save the full table
    df_sorted.to_csv(OUTPUT_PATH, index=False)

    # 5. BEAUTIFUL TERMINAL OUTPUT
    print("\n" + "="*95)
    print(f"{'GENE':<12} | {'MEAN %':<8} | {'DIFF':<6} | {'CPTAC %':<10} | {'AURORA %':<10} | {'CPTAC N':<8} | {'AURORA N':<8}")
    print("-" * 95)

    # Display Top 25 most consistent genes
    for _, row in df_sorted.head(25).iterrows():
        print(f"{row['gene']:<12} | "
              f"{row['Mean_Pct']:<8.1f} | "
              f"{row['Diff']:<6.1f} | "
              f"{row['CPTAC_Percent']:<10.1f} | "
              f"{row['AURORA_Percent']:<10.1f} | "
              f"{int(row['CPTAC_Mutations']):<8} | "
              f"{int(row['AURORA_Mutations']):<8}")
    
    print("-" * 95)
    print(f"📊 Total genes analyzed: {len(df_sorted)}")
    print(f"✅ Final scientific report saved to: {OUTPUT_PATH}")
    print("="*95 + "\n")

if __name__ == "__main__":
    generate_clean_report()
