import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def main():
    print("🚀 Запуск финального доказательства теории Convergence...")
    
    # 1. Загружаем данные по Силенсингу из Пункта 1
    # Мы используем файл, который создали при анализе HLA-экспрессии, там уже есть группы
    silencing_file = 'results/hla_comparison_data.csv'
    try:
        silencing_df = pd.read_csv(silencing_file)
    except FileNotFoundError:
        print("Ошибка: Сначала запусти скрипт analyze_hla_downregulation.py!")
        return

    # Нам нужны только уникальные пациенты и их группы (High/Low Silencing)
    patient_groups = silencing_df[['Tumor_Sample_Barcode', 'Group']].drop_duplicates()

    # 2. Загружаем данные CNA по B2M
    cna_file = 'data/data_cna.txt'
    cna = pd.read_csv(cna_file, sep='\t')
    
    # Вытаскиваем только B2M
    b2m_cna = cna[cna['Hugo_Symbol'] == 'B2M'].melt(
        id_vars=['Hugo_Symbol'], var_name='Tumor_Sample_Barcode', value_name='CNA_Value'
    )
    
    # Статус делеции: True если <= -1
    b2m_cna['Is_B2M_Deleted'] = b2m_cna['CNA_Value'] <= -1.0

    # 3. Мерджим данные
    final_df = pd.merge(patient_groups, b2m_cna[['Tumor_Sample_Barcode', 'Is_B2M_Deleted']], on='Tumor_Sample_Barcode')

    # 4. Строим таблицу для теста Фишера
    # Группа (High/Low) vs Удаление (True/False)
    table = pd.crosstab(final_df['Group'], final_df['Is_B2M_Deleted'])
    
    print("\n--- ТАБЛИЦА СОПРЯЖЕННОСТИ ---")
    print(table)
    print("----------------------------")

    # 5. Считаем статистику
    # Тест Фишера идеален для малых выборок и таблиц 2x2
    odds_ratio, p_value = fisher_exact(table)

    # 6. Итоговый отчет
    print(f"\nРЕЗУЛЬТАТЫ АНАЛИЗА:")
    
    for group in table.index:
        total = table.loc[group].sum()
        deleted = table.loc[group, True]
        perc = (deleted / total) * 100
        print(f"Группа {group:15}: {perc:5.1f}% имеют делецию B2M ({deleted} из {total})")

    print(f"\nP-value (Fisher's exact test): {p_value:.4f}")

    if p_value < 0.05:
        print("\n✅ ГИПОТЕЗА ПОДТВЕРЖДЕНА! Рак выбирает геномное удаление, если не смог спрятать мутации.")
    else:
        print("\n⚠️ Тренд есть, но P-value > 0.05. Возможно, выборка в 122 пациента маловата для этого датасета.")

    # Сохраняем для отчета
    final_df.to_csv('results/convergence_final_data.csv', index=False)

if __name__ == "__main__":
    main()
