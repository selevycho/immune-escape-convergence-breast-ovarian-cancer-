import pandas as pd
import matplotlib.pyplot as plt

def main():
    print("Loading data for pie chart...")
    # Читаем данные, которые сгенерировал первый скрипт
    df = pd.read_csv('results/silencing_summary.csv')

    # Настройки графика
    labels = df['Status']
    sizes = df['Count']
    colors = ['#E64B35', '#4DBBD5'] # Красный для Silenced, Синий для Expressed
    explode = (0.05, 0) # Слегка выдвигаем кусок "Silenced" для акцента

    plt.figure(figsize=(8, 8))
    
    # Отрисовка
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, 
            autopct='%1.1f%%', shadow=False, startangle=140,
            textprops={'fontsize': 14, 'fontweight': 'bold'})
    
    plt.title('Transcriptomic Mutation Silencing in Breast Cancer', 
              pad=20, fontsize=16, fontweight='bold')

    # Сохранение в высоком разрешении
    output_path = 'results/silencing_pie_chart.png'
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, transparent=False, facecolor='white')
    print(f"Pie chart successfully saved to: {output_path}")

if __name__ == "__main__":
    main()
