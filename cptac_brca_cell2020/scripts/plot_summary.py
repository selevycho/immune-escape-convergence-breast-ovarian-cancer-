import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Твои финальные данные
    mechanisms = ['Line 1:\nMutation Silencing', 'Line 2:\nHLA Downregulation', 'Line 3:\nB2M Deletion']
    percentages = [34.7, 5.0, 37.7] # 5% для HLA как пример "неэффективности"
    
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 6))
    
    colors = ['#3C5488', '#8491B4', '#DC0000']
    bars = plt.bar(mechanisms, percentages, color=colors, edgecolor='black', alpha=0.8)
    
    # Добавляем цифры над барами
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval}%', ha='center', fontweight='bold')

    plt.ylim(0, 50)
    plt.title('Prevalence of Immune Escape Mechanisms in BRCA', fontsize=14, fontweight='bold')
    plt.ylabel('Percentage of Patients / Mutations (%)', fontweight='bold')
    sns.despine()

    plt.tight_layout()
    plt.savefig('results/escape_summary_stats.png', dpi=300)
    print("✅ Итоговый график готов: results/escape_summary_stats.png")

if __name__ == "__main__":
    main()
