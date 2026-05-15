import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
DATA_FILE = os.path.join(HOME, "results/line3_deletion_summary.csv")
PLOT_FILE = os.path.join(HOME, "results/line3_aurora_convergence.png")

def main():
    if not os.path.exists(DATA_FILE):
        print("[-] Файл данных не найден. Запусти сначала скрипт статистики.")
        return

    df = pd.read_csv(DATA_FILE)
    df['Del_Freq'] *= 100 

    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    
    palette = {"Low Silencing": "#e74c3c", "High Silencing": "#3498db"}
    ax = sns.barplot(x='Gene', y='Del_Freq', hue='Group', data=df, palette=palette)

    plt.title('Line 3: Genomic Deletions Frequency (AURORA Dataset)', fontsize=14, pad=15)
    plt.ylabel('Deletion Frequency (%)')
    plt.ylim(0, 100)
    
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.1f}%', 
                   (p.get_x() + p.get_width() / 2., p.get_height()), 
                   ha='center', va='center', xytext=(0, 9), textcoords='offset points',
                   fontsize=10, fontweight='bold')

    plt.savefig(PLOT_FILE, dpi=300, bbox_inches='tight')
    print(f"[OK] График сохранен: {PLOT_FILE}")

if __name__ == "__main__":
    main()
