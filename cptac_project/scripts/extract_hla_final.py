import pandas as pd
import os

def main():
    # Загружаем промежуточную матрицу v3, которую мы уже создали
    df = pd.read_csv('results/hla_matrix_v3.tsv', sep='\t')
    
    # Оставляем только нужные проекты
    target_projects = ['Breast_Cancer', 'Ovarian_Cancer']
    df_filtered = df[df['Project'].isin(target_projects)].copy()
    
    # Сохраняем чистую матрицу
    df_filtered.to_csv('results/hla_matrix_final.tsv', sep='\t', index=False)
    
    print("--- Final Matrix Filtered ---")
    print(df_filtered['Project'].value_counts())
    print("Saved to results/hla_matrix_final.tsv")

if __name__ == "__main__":
    main()
