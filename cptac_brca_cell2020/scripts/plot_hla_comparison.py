import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    df = pd.read_csv('results/hla_comparison_data.csv')
    
    # Строгий научный стиль
    sns.set_theme(style="ticks", context="talk")
    plt.figure(figsize=(10, 6.5))

    # Рисуем график
    ax = sns.barplot(
        data=df, 
        x='Hugo_Symbol', 
        y='Expression_FPKM', 
        hue='Group',
        palette={'Low Silencing': '#E64B35', 'High Silencing': '#4DBBD5'}, # Красный и синий
        capsize=.1, 
        err_kws={'linewidth': 1.5},
        edgecolor=".2" # Добавляет тонкую рамку вокруг колонок для четкости
    )

    # Настройка текста и осей
    plt.title('MHC-I Complex Expression by Mutation Silencing Strategy', fontsize=16, fontweight='bold', pad=20)
    plt.ylabel('Expression (FPKM)', fontweight='bold', fontsize=14)
    plt.xlabel('Gene', fontweight='bold', fontsize=14)
    
    # Настройка легенды
    plt.legend(title='Patient Strategy', title_fontsize='13', fontsize='12', frameon=True)
    
    # Убираем верхнюю и правую рамки
    sns.despine()

    output = 'results/hla_downregulation_plot.png'
    plt.tight_layout()
    plt.savefig(output, dpi=300, transparent=False, facecolor='white')
    print(f"\n✅ График успешно сохранен в: {output}")

if __name__ == "__main__":
    main()
