import pandas as pd
import os

# 1. Paths
HOME = os.path.expanduser("~")
CPTAC_PATH = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")
AURORA_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")
OUTPUT_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/final_high_silencing_report_v2.csv")

def generate_verified_report():
    if not os.path.exists(CPTAC_PATH) or not os.path.exists(AURORA_PATH):
        print("[-] Error: Missing input files.")
        return

    # 2. Load CPTAC - берем все гены из этого файла, так как он уже отфильтрован по смыслу
    cptac = pd.read_csv(CPTAC_PATH)
    # Используем >= 50, чтобы не терять пограничные значения
    cptac_high = cptac[cptac['Silenced_Percent'] >= 50].copy()
    cptac_high = cptac_high[['Hugo_Symbol', 'Silenced_Percent', 'Total_Mutations']]
    cptac_high.columns = ['gene', 'CPTAC_Percent', 'CPTAC_N']

    # 3. Load AURORA
    aurora = pd.read_csv(AURORA_PATH)
    aurora_high = aurora[['gene', 'ratio', 'N_mutations']].copy()
    aurora_high.columns = ['gene', 'AURORA_Percent', 'AURORA_N']

    # 4. Merge (Intersection)
    combined = pd.merge(cptac_high, aurora_high, on='gene', how='inner')

    # 5. Metrics
    combined['Mean_Pct'] = (combined['CPTAC_Percent'] + combined['AURORA_Percent']) / 2
    combined['Diff'] = abs(combined['CPTAC_Percent'] - combined['AURORA_Percent'])

    # 6. Sorting
    combined = combined.sort_values(by='Mean_Pct', ascending=False)
    combined.to_csv(OUTPUT_PATH, index=False)

    # 7. BEAUTIFUL TERMINAL OUTPUT
    print("\n" + "="*95)
    print(f"{'GENE':<12} | {'MEAN %':<8} | {'DIFF':<6} | {'CPTAC %':<10} | {'AURORA %':<10} | {'CPTAC N':<8} | {'AURORA N':<8}")
    print("-" * 95)

    for _, row in combined.iterrows():
        print(f"{row['gene']:<12} | "
              f"{row['Mean_Pct']:<8.1f} | "
              f"{row['Diff']:<6.1f} | "
              f"{row['CPTAC_Percent']:<10.1f} | "
              f"{row['AURORA_Percent']:<10.1f} | "
              f"{int(row['CPTAC_N']):<8} | "
              f"{int(row['AURORA_N']):<8}")
    
    print("-" * 95)
    print(f"📊 Total shared high-silencing genes: {len(combined)}")
    print(f"✅ Verified report saved to: {OUTPUT_PATH}")
    print("="*95 + "\n")

if __name__ == "__main__":
    generate_verified_report()
