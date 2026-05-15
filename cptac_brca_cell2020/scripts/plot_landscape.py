import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def main():
    # 1. Собираем данные (пересчитываем для точности)
    mut = pd.read_csv('data/data_mutations.txt', sep='\t', low_memory=False)
    exp = pd.read_csv('data/data_mrna_seq_fpkm.txt', sep='\t')
    cna = pd.read_csv('data/data_cna.txt', sep='\t')

    # Очистка
    valid_muts = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site']
    mut = mut[mut['Variant_Classification'].isin(valid_muts)].drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])
    if 'Entrez_Gene_Id' in exp.columns: exp = exp.drop(columns=['Entrez_Gene_Id'])
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')

    # Расчет Ratio
    merged = pd.merge(mut[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_melted, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    merged['Is_Silenced'] = np.where(merged['Expression_FPKM'] < 1.0, 1, 0)
    df = merged.groupby('Tumor_Sample_Barcode')['Is_Silenced'].mean().reset_index()
    df.columns = ['Tumor_Sample_Barcode', 'Silencing_Ratio']

    # Добавляем B2M статус
    b2m = cna[cna['Hugo_Symbol'] == 'B2M'].melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='CNA')
    b2m['B2M_Status'] = np.where(b2m['CNA'] <= -1, 'Deleted', 'Intact')
    
    final = pd.merge(df, b2m[['Tumor_Sample_Barcode', 'B2M_Status']], on='Tumor_Sample_Barcode')
    final = final.sort_values('Silencing_Ratio', ascending=False) # Сортируем для красоты

    # 2. ВИЗУАЛИЗАЦИЯ
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True, 
                                   gridspec_kw={'height_ratios': [4, 1]})

    # Верхний график: Silencing Ratio
    colors = ['#4DBBD5' if s == 'Intact' else '#E64B35' for s in final['B2M_Status']]
    ax1.bar(range(len(final)), final['Silencing_Ratio'], color=colors, width=0.8)
    ax1.set_ylabel('Mutation Silencing Ratio', fontweight='bold')
    ax1.set_title('BRCA Immune Escape Landscape (N=122)', fontsize=16, fontweight='bold')
    ax1.axhline(final['Silencing_Ratio'].median(), color='black', linestyle='--', alpha=0.5, label='Median Ratio')

    # Нижний график: B2M Status Heatmap
    status_map = {'Deleted': 1, 'Intact': 0}
    status_colors = final['B2M_Status'].map(status_map).values.reshape(1, -1)
    ax2.imshow(status_colors, cmap='Reds', aspect='auto', interpolation='nearest')
    ax2.set_yticks([])
    ax2.set_xlabel('Individual Patients', fontweight='bold')
    ax2.set_ylabel('B2M Status', fontweight='bold')

    # Легенда
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#E64B35', label='B2M Deleted'),
                       Patch(facecolor='#4DBBD5', label='B2M Intact')]
    ax1.legend(handles=legend_elements)

    plt.tight_layout()
    plt.savefig('results/immune_escape_landscape.png', dpi=300)
    print("✅ Ландшафт готов: results/immune_escape_landscape.png")

if __name__ == "__main__":
    main()
