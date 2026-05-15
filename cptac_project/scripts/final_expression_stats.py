import pandas as pd
from scipy import stats
import numpy as np

def main():
    # 1. Load the unified HLA matrix
    df = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    
    # 2. Filter for Target Groups ONLY
    targets = ['Breast_Cancer', 'Ovarian_Cancer']
    df = df[df['Project'].isin(targets)].copy()
    
    genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    results = []

    print(f"Analyzing expression for {len(df)} targeted samples...")

    for gene in genes:
        brca_vals = df[df['Project'] == 'Breast_Cancer'][gene]
        ov_vals = df[df['Project'] == 'Ovarian_Cancer'][gene]
        
        # Calculate Means
        mean_brca = brca_vals.mean()
        mean_ov = ov_vals.mean()
        
        # T-test for significance
        t_stat, p_val = stats.ttest_ind(brca_vals, ov_vals, equal_var=False)
        
        # Log2 Fold Change (OV / BRCA)
        l2fc = np.log2(mean_ov / mean_brca)
        
        results.append({
            'Gene': gene,
            'Mean_BRCA': round(mean_brca, 2),
            'Mean_OV': round(mean_ov, 2),
            'Log2FC': round(l2fc, 3),
            'P_Value': round(p_val, 5),
            'Significant': 'Yes' if p_val < 0.05 else 'No'
        })

    # Save Stats Table
    stats_df = pd.DataFrame(results)
    stats_df.to_csv('results/final_dea_stats_brca_ov.tsv', sep='\t', index=False)
    
    print("\n--- FINAL DEA RESULTS (OV vs BRCA) ---")
    print(stats_df)
    print("\nStats saved to results/final_dea_stats_brca_ov.tsv")

if __name__ == "__main__":
    main()
