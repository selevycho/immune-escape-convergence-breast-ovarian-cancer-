import pandas as pd
import numpy as np

def main():
    print("Loading data...")
    mut_file = 'data/data_mutations.txt'
    exp_file = 'data/data_mrna_seq_fpkm.txt'

    # 1. Read mutations
    mut = pd.read_csv(mut_file, sep='\t', low_memory=False)
    mut = mut[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification', 'HGVSp_Short']]

    # Filter: keep only protein-altering mutations
    valid_mutations = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']
    mut = mut[mut['Variant_Classification'].isin(valid_mutations)]

    # 2. Read RNA expression
    exp = pd.read_csv(exp_file, sep='\t')
    if 'Entrez_Gene_Id' in exp.columns:
        exp = exp.drop(columns=['Entrez_Gene_Id'])

    print("Processing expression matrix...")
    # Melt the matrix to match patients and genes
    exp_melted = exp.melt(id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='Expression_FPKM')

    # 3. Merge mutations with their RNA expression
    print("Identifying expressed neoantigens (Remain vs Disappear)...")
    neo_df = pd.merge(mut, exp_melted, on=['Hugo_Symbol', 'Tumor_Sample_Barcode'], how='inner')

    # 4. Analysis: Determine status based on expression threshold
    # FPKM >= 1.0 means the mutated gene is actively transcribed
    neo_df['Status'] = np.where(neo_df['Expression_FPKM'] >= 1.0, 'Remain', 'Disappear')
    
    # Sort by expression level and save
    neo_df = neo_df.sort_values(by='Expression_FPKM', ascending=False)
    neo_df.to_csv('results/neoantigen_candidates.tsv', sep='\t', index=False)

    # 5. Output Statistics
    total = len(neo_df)
    remain = len(neo_df[neo_df['Status'] == 'Remain'])
    disappear = len(neo_df[neo_df['Status'] == 'Disappear'])

    print("\n--- NEOANTIGEN STATISTICS ---")
    print(f"Total potential neoantigens (mutations) found: {total}")
    if total > 0:
        print(f"EXPRESSED (Remain): {remain} ({round(remain/total*100, 1)}%)")
        print(f"SILENCED (Disappear): {disappear} ({round(disappear/total*100, 1)}%)")
    
    print("\nTop 5 ideal candidates for NetMHCpan (highest FPKM expression):")
    print(neo_df.head(5)[['Hugo_Symbol', 'HGVSp_Short', 'Tumor_Sample_Barcode', 'Expression_FPKM']])

if __name__ == "__main__":
    main()
