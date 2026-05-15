import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
DATA_FILE = os.path.join(HOME, "results/aurora_patient_landscape_data.csv")

def create_plot(df, gene_name):
    # Сортировка по Ratio
    df_sorted = df.sort_values('ratio', ascending=False).reset_index(drop=True)
    n = len(df_sorted)
    status_col = f'{gene_name}_status'

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 9), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
    plt.subplots_adjust(hspace=0.02)

    # 1. Верхний график (Bars)
    colors = ['#e74c3c' if d else '#4eb8ce' for d in df_sorted[status_col]]
    ax1.bar(range(n), df_sorted['ratio'], color=colors, width=0.85)
    ax1.axhline(df_sorted['ratio'].mean(), color='gray', linestyle='--', alpha=0.6)
    
    ax1.set_title(f'AURORA BRCA Immune Escape: {gene_name} Status (N={n})', fontsize=20, fontweight='bold')
    ax1.set_ylabel('Mutation Silencing Ratio', fontsize=14)
    
    # 2. Нижний график (Heatmap)
    status_map = df_sorted[status_col].values.reshape(1, -1)
    sns.heatmap(status_map, cmap='Reds', cbar=False, ax=ax2, linewidths=0.1, linecolor='white')
    
    ax2.set_yticks([])
    ax2.set_ylabel(f'{gene_name}\nStatus', fontsize=12, rotation=0, labelpad=40)
    ax2.set_xlabel('Individual Patients (Sorted by Silencing)', fontsize=14)

    output_name = os.path.join(HOME, f"results/landscape_aurora_{gene_name}.png")
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[OK] График для {gene_name} сохранен.")

def main():
    if not os.path.exists(DATA_FILE):
        print("[-] Сначала запусти скрипт обработки данных!")
        return
    
    df = pd.read_csv(DATA_FILE)
    for gene in ['B2M', 'HLA_A', 'HLA_C']:
        create_plot(df, gene)

if __name__ == "__main__":
    main()
