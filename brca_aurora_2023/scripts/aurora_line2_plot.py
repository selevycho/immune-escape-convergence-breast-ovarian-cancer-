import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
DATA_FILE = os.path.join(HOME, "results/line2_comparison_data.csv")
PLOT_FILE = os.path.join(HOME, "results/line2_hla_downregulation.png")

def plot_line2():
    if not os.path.exists(DATA_FILE):
        print("[-] Error: Comparison data not found. Run stats script first.")
        return

    df = pd.read_csv(DATA_FILE)
    
    # Переводим таблицу в длинный формат для Seaborn
    df_melted = df.melt(id_vars=['Tumor_Sample_Barcode', 'Group'], 
                        value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
                        var_name='Gene', value_name='Z_score')

    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    # Строим график
    ax = sns.violinplot(x='Gene', y='Z_score', hue='Group', data=df_melted, 
                        split=True, inner="quart", palette="Set2")
    
    plt.title('Line 2: HLA/B2M Expression vs Line 1 Silencing Profile (AURORA)', fontsize=14)
    plt.axhline(0, linestyle='--', color='black', alpha=0.3)
    plt.ylabel('Expression Z-score')
    plt.xlabel('Antigen Presentation Machinery')
    
    # Добавляем пояснение
    plt.annotate('Lower values = Stronger Immune Escape', xy=(0.5, -2), 
                 xycoords='data', color='red', fontweight='bold')

    plt.savefig(PLOT_FILE, dpi=300)
    print(f"[+] Plot saved to: {PLOT_FILE}")
    plt.show()

if __name__ == "__main__":
    plot_line2()
