import pandas as pd
import matplotlib.pyplot as plt

def main():
    # 1. Читаем файл с мутациями
    try:
        df = pd.read_csv('results/detailed_mutations.tsv', sep='\t')
    except FileNotFoundError:
        print("Error: results/detailed_mutations.tsv not found.")
        return

    # Переименовываем колонки для красивого отображения в дипломе
    df.columns = ['Patient ID', 'Cancer Type', 'Mutated Gene', 'Mutation Type', 'Protein Change']

    # 2. Настраиваем фигуру
    fig, ax = plt.subplots(figsize=(10, 2.5)) # Делаем её невысокой, так как строк всего 4
    ax.axis('off')

    # 3. Рисуем таблицу
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        colColours=["#e0e0e0"] * len(df.columns) # Серый фон для заголовков
    )

    # 4. Форматирование текста и ячеек
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 2)

    # Добавляем цвета строк в зависимости от типа рака
    for i in range(len(df)):
        cancer_type = df.iloc[i]['Cancer Type']
        # Розовый для груди, Голубой для яичников
        bg_color = "#ffe6e6" if cancer_type == 'Breast_Cancer' else "#e6f3ff"
        for j in range(len(df.columns)):
            table[(i + 1, j)].set_facecolor(bg_color)

    plt.title('Identified Somatic Mutations in HLA/B2M Complex', 
              fontsize=14, fontweight='bold', pad=15)

    # 5. Сохранение
    plt.tight_layout()
    plt.savefig('results/mutation_summary_table.png', dpi=300, bbox_inches='tight')
    print("Success! Mutation table saved to results/mutation_summary_table.png")

if __name__ == "__main__":
    main()
