import pandas as pd
import os

# 1. Пути к файлам
HOME = os.path.expanduser("~")
# Это файл, где соединены мутации и Z-scores (то, что я называл Merged Data)
# Если ты его переименовал, проверь название. Обычно это был aurora_mutations_with_zscore.csv
FULL_DATA_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/aurora_mutations_with_zscore.csv")

# Твои два списка-фильтра
HIGH_LIST_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")
LOW_LIST_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/low_silencing_aurora_all.csv")

def create_detailed_table(target_genes_path, output_name):
    if not os.path.exists(FULL_DATA_PATH):
        print(f"[-] Ошибка: Не найден файл с полными данными {FULL_DATA_PATH}")
        return

    # Загружаем полные данные и список нужных генов
    full_df = pd.read_csv(FULL_DATA_PATH)
    target_genes = pd.read_csv(target_genes_path)['gene'].unique()

    # Фильтруем: только нужные гены и только ЗАГЛУШЕННЫЕ (Z < -1.0)
    silenced = full_df[(full_df['Hugo_Symbol'].isin(target_genes)) & (full_df['Z_score'] < -1.0)].copy()

    # Создаем сводную таблицу (Pivot)
    detailed = pd.crosstab(silenced['Hugo_Symbol'], silenced['Variant_Classification'])

    # Список колонок как на твоем скриншоте (чтобы порядок был идеальный)
    cols = [
        'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
        'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 
        'Splice_Site', 'Translation_Start_Site'
    ]

    # Добавляем пустые колонки, если таких типов мутаций не нашлось
    for c in cols:
        if c not in detailed.columns:
            detailed[c] = 0

    # Оставляем только нужные и считаем Total
    detailed = detailed[cols]
    detailed['Total_Silenced'] = detailed.sum(axis=1)
    
    # Сортируем по убыванию общего количества
    detailed = detailed.sort_values(by='Total_Silenced', ascending=False).reset_index()

    # Сохраняем
    out_path = os.path.join(HOME, f"immune_escape_project/brca_aurora_2023/results/{output_name}")
    detailed.to_csv(out_path, index=False)
    print(f"[+] Таблица сохранена: {out_path}")
    print(detailed.head(10).to_string(index=False))
    print("-" * 100)

if __name__ == "__main__":
    # Делаем детализацию для High Silencing
    create_detailed_table(HIGH_LIST_PATH, "detailed_high_silencing_aurora.csv")
    # Делаем детализацию для Low Silencing
    create_detailed_table(LOW_LIST_PATH, "detailed_low_silencing_aurora.csv")
