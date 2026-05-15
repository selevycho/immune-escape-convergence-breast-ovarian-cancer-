import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
DATA_FILE = os.path.join(HOME, "results/aurora_patient_landscape_data.csv")

def create_perfect_plot(df, gene_name):
    # 1. Сортировка
    df_sorted = df.sort_values('ratio', ascending=False).reset_index(drop=True)
    n = len(df_sorted)
    status_col = f'{gene_name}_status'

    # 2. Создаем фигуру
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
    plt.subplots_adjust(hspace=0.05)

    # Цвета: Красный - Deletion (1), Голубой - Intact (0)
    colors = ['#e74c3c' if d == 1 else '#4eb8ce' for d in df_sorted[status_col]]
    
    # 3. Верхний график: Mutation Silencing Ratio
    # Используем align='edge' и фиксированный xlim для идеального выравнивания
    x_pos = np.arange(n)
    ax1.bar(x_pos, df_sorted['ratio'], color=colors, width=0.8, align='center')
    
    # ФИКСИРУЕМ ШКАЛУ Y для всех графиков (от 0 до 0.45)
    ax1.set_ylim(0, 0.45)
    ax1.axhline(df_sorted['ratio'].mean(), color='black', linestyle='--', alpha=0.3, label='Mean Ratio')
    
    ax1.set_title(f'AURORA BRCA Immune Escape Landscape: {gene_name} Analysis (N={n} patients)', fontsize=22, fontweight='bold', pad=20)
    ax1.set_ylabel('Mutation Silencing Ratio\n(Line 1: RNA suppression)', fontsize=14)

    # ДОБАВЛЯЕМ ЛЕГЕНДУ
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label=f'{gene_name} Deleted (Genomic Loss)'),
        Patch(facecolor='#4eb8ce', label=f'{gene_name} Intact (Genomic Presence)')
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=14, frameon=True)

    # 4. Нижний график: Heatmap статуса
    # Чтобы heatmap совпал с барами, используем imshow или pcolormesh
    status_data = df_sorted[status_col].values.reshape(1, -1)
    # Используем кастомную карту цветов: 0 -> белый/голубой, 1 -> темно-красный
    cmap = sns.dark_palette("#e74c3c", as_cmap=True)
    
    sns.heatmap(status_data, cmap=['#fdfefe', '#800000'], cbar=False, ax=ax2, 
                linewidths=0.5, linecolor='#ecf0f1')
    
    ax2.set_yticks([])
    ax2.set_ylabel(f'{gene_name}\nStatus', fontsize=14, rotation=0, labelpad=45, verticalalignment='center')
    ax2.set_xlabel('Individual Patients (Ranked by Silencing Efficiency)', fontsize=14, labelpad=15)
    
    # Убираем лишние цифры с оси X, оставляем только основные деления
    ax2.set_xticks(np.arange(0, n, 5) + 0.5)
    ax2.set_xticklabels(np.arange(0, n, 5))

    # Жестко фиксируем границы осей, чтобы они не разъезжались
    ax1.set_xlim(-0.5, n - 0.5)
    ax2.set_xlim(0, n)

    output_path = os.path.join(HOME, f"results/landscape_aurora_{gene_name}_v3.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[OK] Исправленный график {gene_name} сохранен.")

def main():
    if not os.path.exists(DATA_FILE):
        print("[-] Ошибка: Данные не найдены.")
        return
    df = pd.read_csv(DATA_FILE)
    for gene in ['B2M', 'HLA_A', 'HLA_C']:
        create_perfect_plot(df, gene)

if __name__ == "__main__":
    main()
