import pandas as pd
import numpy as np

def main():
    # 1. Загружаем данные
    expr_df = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    mut_df = pd.read_csv('results/final_brca_ov_mutations.tsv', sep='\t')
    
    # ОЧИСТКА: убираем пробелы и приводим к одному регистру
    expr_df['Case_ID'] = expr_df['Case_ID'].astype(str).str.strip().str.upper()
    mut_df['Case_ID'] = mut_df['Case_ID'].astype(str).str.strip().str.upper()
    
    mutant_cases = mut_df['Case_ID'].unique()
    
    print(f"DEBUG: Unique Case IDs in Expression: {len(expr_df['Case_ID'].unique())}")
    print(f"DEBUG: Unique Case IDs in Mutations: {len(mutant_cases)}")
    print(f"DEBUG: Samples of Mutation IDs: {mutant_cases[:5]}")

    # 2. Помечаем статус
    expr_df['Status'] = expr_df['Case_ID'].apply(lambda x: 'Mutant' if x in mutant_cases else 'Wild_Type')
    
    # Проверка: сколько в итоге совпало?
    mutants_found = expr_df[expr_df['Status'] == 'Mutant']
    print(f"DEBUG: Matches found in Expression table: {len(mutants_found)}")

    if len(mutants_found) == 0:
        print("\nERROR: No matches found again. Let's try partial matching...")
        # Попробуем найти вхождения (если ID в одной таблице длиннее, чем в другой)
        def partial_match(cid):
            for m_id in mutant_cases:
                if m_id in cid or cid in m_id:
                    return 'Mutant'
            return 'Wild_Type'
        expr_df['Status'] = expr_df['Case_ID'].apply(partial_match)
        mutants_found = expr_df[expr_df['Status'] == 'Mutant']
        print(f"DEBUG: Matches found after partial matching: {len(mutants_found)}")

    genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    
    print("\n" + "="*60)
    print("FINAL EXPRESSION COMPARISON: MUTANTS VS WILD TYPE")
    print("="*60)
    
    results = []
    for gene in genes:
        mean_mutant = expr_df[expr_df['Status'] == 'Mutant'][gene].mean()
        mean_wild = expr_df[expr_df['Status'] == 'Wild_Type'][gene].mean()
        ratio = (mean_mutant / mean_wild * 100) if mean_wild > 0 else 0
        
        results.append({
            'Gene': gene,
            'Mean_WildType': round(mean_wild, 2),
            'Mean_Mutant': round(mean_mutant, 2) if not np.isnan(mean_mutant) else "N/A",
            'Retention (%)': round(ratio, 1) if not np.isnan(ratio) else "N/A"
        })
    
    print(pd.DataFrame(results))
    
    if not mutants_found.empty:
        print("\n--- DATA FOR THE MUTANT PATIENTS ---")
        print(mutants_found[['Case_ID', 'Project'] + genes])
    
    expr_df.to_csv('results/mutant_impact_final.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
