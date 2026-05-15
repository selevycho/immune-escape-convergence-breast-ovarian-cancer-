import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    print("Loading data for plot...")
    df = pd.read_csv('results/hla_groups_expression.tsv', sep='\t')

    sns.set_theme(style="whitegrid", context="talk")
    plt.figure(figsize=(10, 6))

    # ИСПРАВЛЕНО: errwidth заменен на err_kws
    ax = sns.barplot(
        data=df, 
        x='Hugo_Symbol', 
        y='Expression_FPKM', 
        hue='Group',
        palette=["#4DBBD5", "#E64B35"], # Синий и Красный
        capsize=.1,
        err_kws={'linewidth': 1.5}
    )

    plt.title('MHC-I Complex Expression by Mutation Silencing Strategy', pad=20, fontweight='bold')
    plt.xlabel('Gene', fontweight='bold')
    plt.ylabel('Expression (FPKM)', fontweight='bold')
    
    plt.legend(title='Mutation Silencing Strategy', title_fontsize='13', fontsize='12')
    sns.despine()

    output_path = 'results/hla_escape_barplot.png'
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Plot successfully saved to: {output_path}")

if __name__ == "__main__":
    main()
