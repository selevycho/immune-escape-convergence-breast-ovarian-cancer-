import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Пути к твоим файлам на сервере
HOME = "/home/fr/fr_fr/fr_os136/immune_escape_project/brca_aurora_2023"
MERGED_DATA = os.path.join(HOME, "results/line1_aurora_merged_data.csv")
CNA_FILE = os.path.join(HOME, "data/data_log2_cna.txt")
PLOT_OUTPUT = os.path.join(HOME, "results/aurora_immune_landscape_final.png")

def main():
    print("[+] Начинаю полную сборку ландшафта для AURORA...")

    # 1. Загрузка мутаций и расчет Silencing Ratio
    # Мы превращаем Z-score в "флажки" (1 - замолчал, 0 - нет)
    df_muts = pd.read_csv(MERGED_DATA)
    df_muts['Is_Silenced'] = np.where(df_muts['Z_score'] < -1.0, 1, 0)
    
    # Считаем Ratio для каждого пациента: (Сумма единичек) / (Общее кол-во мутаций)
    patient_stats = df_muts.groupby('Tumor_Sample_Barcode').agg(
        total_muts=('Is_Silenced', 'count'),
        silenced_count=('Is_Silenced', 'sum')
    ).reset_index()
    patient_stats['silencing_ratio'] = patient_stats['silenced_count'] / patient_stats['total_muts']

    # 2. Загрузка CNA для B2M (Entrez ID: 567)
    cna = pd.read_csv(CNA_FILE, sep='\t', low_memory=False)
    b2m_data = cna[cna['Entrez_Gene_Id'] == 567].drop(columns=['Entrez_Gene_Id'])
    
    # Переворачиваем таблицу CNA, чтобы приклеить её к пациентам
    b2m_status = b2m_data.melt(var_name='Tumor_Sample_Barcode', value_name='log2_cna')
    # Порог делеции -0.3 (все что ниже — считаем потерей гена)
    b2m_status['b2m_deleted'] = np.where(b2m_status['log2_cna'] < -0.3, 1, 0)

    # 3. Объединяем всё в одну таблицу
    final_df = pd.merge(patient_stats, b2m_status[['Tumor_Sample_Barcode', 'b2m_deleted']], on='Tumor_Sample_Barcode')
    
    # Сортируем пациентов от самого высокого Ratio к самому низкому
    final_df = final_df.sort_values('silencing_ratio', ascending=False).reset_index(drop=True)
    n_patients = len(final_df)

    # 4. РИСУЕМ ГРАФИК (Landscape)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 9), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
    plt.subplots_adjust(hspace=0.02)

    # Верхняя часть: Бары Ratio
    # Если B2M удален (1) — красим в красный, если цел (0) — в голубой
    colors = ['#e74c3c' if d else '#4eb8ce' for d in final_df['b2m_deleted']]
    ax1.bar(range(n_patients), final_df['silencing_ratio'], color=colors, width=0.85)
    
    # Линия среднего Ratio по когорте
    ax1.axhline(final_df['silencing_ratio'].median(), color='black', linestyle='--', alpha=0.5)
    
    ax1.set_title(f'AURORA Metastatic BRCA Immune Escape Landscape (N={n_patients})', fontsize=18, fontweight='bold')
    ax1.set_ylabel('Mutation Silencing Ratio\n(Line 1: RNA)', fontsize=12)
    
    # Легенда
    from matplotlib.patches import Patch
    ax1.legend(handles=[Patch(facecolor='#e74c3c', label='B2M Deleted (Line 3)'),
                        Patch(facecolor='#4eb8ce', label='B2M Intact')], loc='upper right')

    # Нижняя часть: Полоска статуса B2M (Heatmap)
    b2m_map = final_df['b2m_deleted'].values.reshape(1, -1)
    sns.heatmap(b2m_map, cmap='Reds', cbar=False, ax=ax2, linewidths=0.1, linecolor='white')
    
    ax2.set_yticks([])
    ax2.set_ylabel('B2M Status', fontsize=12, rotation=0, labelpad=40)
    ax2.set_xlabel('Individual Patients (Sorted by Silencing)', fontsize=12)

    plt.savefig(PLOT_OUTPUT, dpi=300, bbox_inches='tight')
    print(f"[OK] Всё готово! График сохранен как: {PLOT_OUTPUT}")

if __name__ == "__main__":
    main()
