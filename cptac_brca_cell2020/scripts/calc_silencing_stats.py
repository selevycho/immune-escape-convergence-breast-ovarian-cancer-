import pandas as pd
import numpy as np

def main():
    print("Initializing Ultra-Robust Pipeline...\n")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # --- 1. ЗАГРУЗКА И СТРОГАЯ ФИЛЬТРАЦИЯ ДНК ---
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    total_raw = len(mut)
    
    # Расширенный список мутаций, которые РЕАЛЬНО меняют белок
    valid_mutations = [
        'Missense_Mutation', 'Nonsense_Mutation', 
        'Frame_Shift_Del', 'Frame_Shift_Ins', 
        'In_Frame_Del', 'In_Frame_Ins',
        'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation'
    ]
    mut_filtered = mut[mut['Variant_Classification'].isin(valid_mutations)]
    
    # Убираем дубликаты: если у пациента 2 мутации в одном гене, 
    # РНК этого гена мы должны проверить только один раз!
    unique_muts = mut_filtered.drop_duplicates(subset=['Hugo_Symbol', 'Tumor_Sample_Barcode'])
    
    print(f"[DNA] Total raw mutations found: {total_raw}")
    print(f"[DNA] Protein-altering mutations: {len(mut_filtered)}")
    print(f"[DNA] Unique Gene-Patient targets: {len(unique_muts)}\n")

    # --- 2. ЗАГРУЗКА И ОЧИСТКА РНК ---
    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])
    
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')
    
    # Удаляем битые данные (где РНК вообще не измерялась)
    exp_clean = exp_melted.dropna(subset=['Expression_FPKM'])

    # --- 3. ОБЪЕДИНЕНИЕ (Ищем улики на складе) ---
    neo_df = pd.merge(unique_muts[['Hugo_Symbol', 'Tumor_Sample_Barcode']], exp_clean, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')
    print(f"[MERGE] Valid targets with available RNA data: {len(neo_df)}\n")

    # --- 4. КЛАССИФИКАЦИЯ (Immune Escape) ---
    neo_df['Status'] = np.where(neo_df['Expression_FPKM'] < 1.0, 'Silenced (Disappear)', 'Expressed (Remain)')

    # --- 5. ФИНАЛЬНАЯ СТАТИСТИКА ---
    total_evaluated = len(neo_df)
    counts = neo_df['Status'].value_counts()
    
    silenced = counts.get('Silenced (Disappear)', 0)
    expressed = counts.get('Expressed (Remain)', 0)
    
    silenced_pct = (silenced / total_evaluated) * 100
    expressed_pct = (expressed / total_evaluated) * 100

    print("=== FINAL SCIENTIFIC STATISTICS ===")
    print(f"Total Evaluated Neoantigens: {total_evaluated}")
    print(f"Silenced (Escape successful): {silenced} ({silenced_pct:.1f}%)")
    print(f"Expressed (Escape failed): {expressed} ({expressed_pct:.1f}%)")
    print("===================================\n")

    # Сохраняем для графика
    summary_df = pd.DataFrame({
        'Status': ['Silenced (Disappear)', 'Expressed (Remain)'],
        'Count': [silenced, expressed]
    })
    summary_df.to_csv('results/silencing_summary.csv', index=False)
    print("Verification complete. Stats saved for plotting.")

if __name__ == "__main__":
    main()
