import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
DATA_FILE = os.path.join(HOME, "results/hla_quartile_comparison_aurora.csv")
PLOT_FILE = os.path.join(HOME, "results/hla_quartile_plot_aurora.png")

def main():
    if not os.path.exists(DATA_FILE):
        print("[-] Файл данных не найден. Сначала запусти скрипт статистики.")
        return

    df = pd.read_csv(DATA_FILE)

    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    
    # Цвета: Красный для "открытых" (Low), Синий для "скрытых" (High)
    palette = {"Low Silencing": "#e74c3c", "High Silencing": "#3498db"}
    
    ax = sns.boxplot(x='Hugo_Symbol', y='Z_score_HLA', hue='Group', data=df, 
                     palette=palette, showfliers=False, width=0.6)
    
    sns.stripplot(x='Hugo_Symbol', y='Z_score_HLA', hue='Group', data=df, 
                  dodge=True, alpha=0.4, palette=palette, size=4, ax=ax)

    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.title('Line 2: HLA Expression by Mutation Silencing Quartiles (AURORA)', fontsize=14, pad=20)
    plt.ylabel('Expression Z-score')
    plt.xlabel('MHC-I Components')
    
    # Чистим легенду
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], title='Patient Quartile')

    plt.savefig(PLOT_FILE, dpi=300, bbox_inches='tight')
    print(f"[OK] График сохранен: {PLOT_FILE}")

if __name__ == "__main__":
    main()
