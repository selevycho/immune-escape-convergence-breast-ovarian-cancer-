import pandas as pd
import numpy as np
from scipy import stats

def main():
    df = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    results = []

    brca = df[df['Project'] == 'Breast_Cancer']
    ov = df[df['Project'] == 'Ovarian_Cancer']

    for gene in genes:
        mean_brca = brca[gene].mean()
        mean_ov = ov[gene].mean()
        
        # Log2 Fold Change (если < 0, значит в яичниках ген слабее, чем в груди)
        lfc = np.log2(mean_ov / mean_brca)
        
        # P-value (Mann-Whitney тест, он надежнее для таких данных)
        stat, p_val = stats.mannwhitneyu(ov[gene].dropna(), brca[gene].dropna())
        
        results.append({
            'Gene': gene,
            'Avg_Breast': round(mean_brca, 2),
            'Avg_Ovarian': round(mean_ov, 2),
            'Log2FC_OV_vs_BRCA': round(lfc, 3),
            'P_value': f"{p_val:.4e}"
        })

    res_df = pd.DataFrame(results)
    print("\n--- DIFFERENTIAL ANALYSIS: OV vs BRCA ---")
    print(res_df)
    res_df.to_csv('results/hla_diff_results.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
