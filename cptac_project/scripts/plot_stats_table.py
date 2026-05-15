import pandas as pd
import matplotlib.pyplot as plt

def main():
    # 1. Загружаем данные DEA
    try:
        df = pd.read_csv('results/final_dea_stats_brca_ov.tsv', sep='\t')
    except FileNotFoundError:
        print("Error: File results/final_dea_stats_brca_ov.tsv not found. Run DEA script first.")
        return

    # 2. Подготовка фигуры
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis('off') # Убираем оси координат

    # 3. Создание таблицы
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        colColours=["#f2f2f2"] * len(df.columns) # Серый заголовок
    )

    # 4. Стилизация (делаем шрифт побольше и выделяем значимые строки)
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 2) # Масштабируем высоту ячеек

    # Подсвечиваем ячейки со значимыми результатами (p < 0.05)
    for i in range(len(df)):
        if df.iloc[i]['Significant'] == 'Yes':
            for j in range(len(df.columns)):
                table[(i + 1, j)].set_facecolor("#e6f3ff") # Светло-голубой для значимых генов

    plt.title('Differential Expression Statistics: Ovarian vs Breast Cancer', 
              fontsize=14, fontweight='bold', pad=20)

    # 5. Сохранение
    plt.tight_layout()
    plt.savefig('results/final_stats_table_image.png', dpi=300, bbox_inches='tight')
    print("Stats table image saved to results/final_stats_table_image.png")

if __name__ == "__main__":
    main()
