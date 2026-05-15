import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # Load our matrix
    df = pd.read_csv('results/hla_expression_matrix.tsv', sep='\t')
    
    # We need to reshape data for plotting (Long format)
    # This turns genes from columns into rows
    plot_df = pd.melt(df, id_vars=['Tissue_Type'], value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
                      var_name='Gene', value_name='Expression_TPM')

    # Create the plot
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=plot_df, x='Gene', y='Expression_TPM', hue='Tissue_Type')
    
    plt.title('Comparison of HLA Class I Expression: Tumor vs Normal')
    plt.yscale('log') # Use log scale because B2M is often much higher than HLA
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Save the plot
    output_img = 'results/hla_comparison_plot.png'
    plt.savefig(output_img)
    print(f"Plot saved successfully to: {output_img}")

if __name__ == "__main__":
    main()
