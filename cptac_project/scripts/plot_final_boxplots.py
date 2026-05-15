import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # 1. Load and filter data
    df = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    df = df[df['Project'].isin(['Breast_Cancer', 'Ovarian_Cancer'])]
    
    # 2. Reshape for Seaborn (Melt)
    # This turns the table from "Wide" to "Long" format for easier plotting
    plot_df = df.melt(id_vars=['Case_ID', 'Project'], 
                      value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
                      var_name='Gene', value_name='Expression_TPM')

    # 3. Create the Plot
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    
    ax = sns.boxplot(x='Gene', y='Expression_TPM', hue='Project', 
                     data=plot_df, palette={'Breast_Cancer': '#ff9999', 'Ovarian_Cancer': '#66b3ff'},
                     showfliers=False) # Hide extreme outliers for clarity
    
    # Add individual points (stripplot) to show the distribution
    sns.stripplot(x='Gene', y='Expression_TPM', hue='Project', 
                  data=plot_df, dodge=True, alpha=0.3, jitter=True, 
                  palette={'Breast_Cancer': 'black', 'Ovarian_Cancer': 'black'},
                  legend=False)

    plt.title('Comparison of HLA/B2M Expression: Breast vs Ovarian Cancer', fontsize=15, fontweight='bold')
    plt.ylabel('Expression (TPM)', fontsize=12)
    plt.xlabel('Gene Symbol', fontsize=12)
    
    # Save the Plot
    plt.savefig('results/final_expression_comparison_clean.png', dpi=300)
    print("Final Boxplot saved to results/final_expression_comparison_clean.png")

if __name__ == "__main__":
    main()
