import pandas as pd
import os

# Пути к файлам
HOME = os.path.expanduser("~")
CPTAC_PATH = os.path.join(HOME, "immune_escape_project/cptac_brca_cell2020/results/high_silencing_all.csv")
AURORA_PATH = os.path.join(HOME, "immune_escape_project/brca_aurora_2023/results/high_silencing_aurora_all.csv")

def verify_raw_intersection():
    if not os.path.exists(CPTAC_PATH) or not os.path.exists(AURORA_PATH):
        print("[-] Ошибка: Файлы не найдены.")
        return

    # 1. Загружаем CPTAC
    cptac_df = pd.read_csv(CPTAC_PATH)
    # В CPTAC колонка называется 'Hugo_Symbol'
    cptac_genes = set(cptac_df['Hugo_Symbol'].unique())

    # 2. Загружаем AURORA
    aurora_df = pd.read_csv(AURORA_PATH)
    # В AURORA колонка называется 'gene'
    aurora_genes = set(aurora_df['gene'].unique())

    # 3. Считаем пересечение
    shared_genes = cptac_genes.intersection(aurora_genes)

    print("\n" + "="*50)
    print("RAW GENE LIST VERIFICATION (HIGH SILENCING)")
    print("="*50)
    print(f"Всего генов в CPTAC (Primary):    {len(cptac_genes)}")
    print(f"Всего генов в AURORA (Metastatic): {len(aurora_genes)}")
    print("-" * 50)
    print(f"ОБЩИХ ГЕНОВ (INTERSECTION):        {len(shared_genes)}")
    print("-" * 50)

    if len(shared_genes) > 0:
        print(f"Список общих генов: {sorted(list(shared_genes))}")
    else:
        print("Общих генов не найдено.")
    print("="*50 + "\n")

if __name__ == "__main__":
    verify_raw_intersection()
