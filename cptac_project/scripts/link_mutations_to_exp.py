import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    if not os.path.exists('results/detailed_mutations.tsv'):
        print("Error: detailed_mutations.tsv not found.")
        return
        
    exp = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    mut = pd.read_csv('results/detailed_mutations.tsv', sep='\t')
    
    # Помечаем мутантов
    exp['Status'] = 'Wild_Type'
    mutant_cases = mut['Case_ID'].unique()
    exp.loc[exp['Case_ID'].isin(mutant_cases), 'Status'] = 'Mutant'
    
    # Генерируем график
    plt.figure(figsize=(10, 7))
    sns.set_style("white")

    # Рисуем сначала всех "здоровых" (серым цветом, чтобы не отвлекали)
    sns.stripplot(data=exp[exp['Status']=='Wild_Type'], x='Status', y='B2M', 
                  color='lightgray', alpha=0.5, jitter=0.2, label='Wild Type')

    # Рисуем мутантов ПОВЕРХ (ярко-красным и крупно)
    sns.stripplot(data=exp[exp['Status']=='Mutant'], x='Status', y='B2M', 
                  hue='Project', palette='Set1', size=12, linewidth=2, 
                  edgecolor='black', jitter=0.1)

    plt.axhline(y=exp[exp['Status'] == 'Wild_Type']['B2M'].mean(), 
                color='red', linestyle='--', label='Mean Expression')

    plt.title('B2M Expression: Wild Type vs Specific Mutants', fontsize=14, fontweight='bold')
    plt.ylabel('Expression (TPM)')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.savefig('results/mutation_impact_visual_v2.png', dpi=300)
    print("\nNew chart saved: results/mutation_impact_visual_v2.png")

if __name__ == "__main__":
    main()
