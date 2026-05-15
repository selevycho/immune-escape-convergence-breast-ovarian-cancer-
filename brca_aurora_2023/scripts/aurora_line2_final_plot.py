import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
IN_FILE = os.path.join(HOME, "results/hla_quartile_comparison_final.csv")
OUT_PLOT = os.path.join(HOME, "results/line2_final_plot.png")

def main():
    df = pd.read_csv(IN_FILE)
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    
    palette = {"Low Silencing": "#e74c3c", "High Silencing": "#3498db"}
    
    # Боксплоты + точки
    ax = sns.boxplot(data=df, x='Hugo_Symbol', y='Z_Exp', hue='Group', 
                     palette=palette, showfliers=False)
    sns.stripplot(data=df, x='Hugo_Symbol', y='Z_Exp', hue='Group', 
                  dodge=True, alpha=0.3, palette=palette, ax=ax)
    
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.title('Line 2: HLA/B2M Expression by Silencing Groups (AURORA)')
    plt.ylabel('Expression Z-score')
    
    # Фикс легенды
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], title='Group')
    
    plt.savefig(OUT_PLOT, dpi=300, bbox_inches='tight')
    print(f"[OK] График готов: {OUT_PLOT}")

if __name__ == "__main__":
    main()
