import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    exp = pd.read_csv('results/hla_matrix_final.tsv', sep='\t')
    mut = pd.read_csv('results/detailed_mutations.tsv', sep='\t')
    
    # Делаем таблицу "длинной" для удобства сравнения
    melted = exp.melt(id_vars=['Case_ID', 'Project'], 
                     value_vars=['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
                     var_name='Gene', value_name='TPM')
    
    # Помечаем конкретные мутированные гены у конкретных людей
    melted['Is_Mutant'] = False
    for _, row in mut.iterrows():
        mask = (melted['Case_ID'] == row['Case_ID']) & (melted['Gene'] == row['Gene'])
        melted.loc[mask, 'Is_Mutant'] = True

    # Визуализация
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    # Рисуем все точки (фоновые)
    sns.stripplot(data=melted[melted['Is_Mutant']==False], x='Gene', y='TPM', 
                  color='lightgray', alpha=0.3, jitter=True)
    
    # Рисуем мутантов (ЯРКО)
    sns.stripplot(data=melted[melted['Is_Mutant']==True], x='Gene', y='TPM', 
                  hue='Case_ID', size=15, marker='X', linewidth=2, edgecolor='black')

    plt.title('Direct Impact of HLA/B2M Mutations on Expression', fontsize=14, fontweight='bold')
    plt.ylabel('Expression (TPM)')
    
    plt.savefig('results/hla_specific_mutation_impact.png', dpi=300)
    print("Chart saved: results/hla_specific_mutation_impact.png")
    
    # Вывод в консоль для проверки
    print("\nCheck these values:")
    print(melted[melted['Is_Mutant']==True][['Case_ID', 'Gene', 'TPM']])

if __name__ == "__main__":
    main()
