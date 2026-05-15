import pandas as pd
import os

# 1. Пути к твоим файлам
HOME = os.path.expanduser("~")
MERGED_DATA = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/line1_aurora_merged_data.csv")
HIGH_LIST = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")
LOW_LIST = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/low_silencing_aurora_all.csv")
RESULT_DIR = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results")

def create_detailed_table(genes_filter, output_name):
    print(f"\n[+] Формирование таблицы: {output_name}")
    
    # Загружаем данные
    df = pd.read_csv(MERGED_DATA)
    df.columns = [c.strip() for c in df.columns]

    # --- АВТОМАТИЧЕСКИЙ ПОИСК КОЛОНОК ---
    # Ищем колонку с геном (Hugo_Symbol, gene, и т.д.)
    gene_col = next((c for c in df.columns if any(x in c.lower() for x in ['hugo', 'symbol', 'gene'])), None)
    
    # Ищем колонку с типом мутации (ищем ту, где есть слова Missense или Nonsense)
    var_col = None
    for col in df.columns:
        unique_vals = str(df[col].unique())
        if 'Mutation' in unique_vals or 'Sense' in unique_vals or 'Shift' in unique_vals:
            var_col = col
            break
    
    if not gene_col or not var_col:
        print(f"[-] ОШИБКА: Не удалось распознать колонки. Доступные: {list(df.columns)}")
        return

    print(f"[i] Найдено: Гены -> '{gene_col}', Тип мутации -> '{var_col}'")

    # --- ФИЛЬТРАЦИЯ И ПОДСЧЕТ ---
    # Оставляем только нужные гены и только заглушенные (Z < -1.0)
    subset = df[(df[gene_col].isin(genes_filter)) & (df['Z_score'] < -1.0)].copy()
    
    if subset.empty:
        print(f"[!] Нет данных для заглушенных мутаций в этих генах.")
        return

    # Создаем сводную таблицу (Pivot)
    pivot = pd.crosstab(subset[gene_col], subset[var_col])

    # Список колонок из твоего ТЗ (скриншота)
    target_cols = [
        'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
        'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 
        'Splice_Site', 'Translation_Start_Site'
    ]

    # Гарантируем наличие всех колонок
    for col in target_cols:
        if col not in pivot.columns:
            pivot[col] = 0

    # Оставляем только нужные, считаем Total и сортируем
    pivot = pivot[target_cols]
    pivot['Total_Silenced'] = pivot.sum(axis=1)
    pivot = pivot.sort_values(by='Total_Silenced', ascending=False).reset_index()
    pivot = pivot.rename(columns={gene_col: 'Hugo_Symbol'})

    # Сохраняем результат
    out_path = os.path.join(RESULT_DIR, output_name)
    pivot.to_csv(out_path, index=False)
    print(f"[OK] Сохранено в: {out_path}")
    print(pivot.head(10).to_string(index=False))

if __name__ == "__main__":
    if os.path.exists(MERGED_DATA):
        # Загружаем списки генов
        high_genes = pd.read_csv(HIGH_LIST)['gene'].unique()
        low_genes = pd.read_csv(LOW_LIST)['gene'].unique()
        
        # Создаем детальные отчеты
        create_detailed_table(high_genes, "detailed_silenced_high_aurora.csv")
        create_detailed_table(low_genes, "detailed_silenced_low_aurora.csv")
    else:
        print(f"[-] Файл {MERGED_DATA} не найден. Проверь путь.")
