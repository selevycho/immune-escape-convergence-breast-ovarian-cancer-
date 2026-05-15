import pandas as pd
import gzip
import os

def main():
    # 1. Загружаем инвентарь
    if not os.path.exists('results/maf_file_inventory.tsv'):
        print("Error: Inventory file not found!")
        return
        
    inv = pd.read_csv('results/maf_file_inventory.tsv', sep='\t')
    inv = inv[inv['Cancer_Type'].isin(['Breast_Cancer', 'Ovarian_Cancer'])]
    
    target_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    mutations = []

    print(f"Checking folders for {len(inv)} cases...")

    for _, row in inv.iterrows():
        folder = row['Folder_ID']
        case = row['Case_ID']
        cancer = row['Cancer_Type']
        
        # ПРОВЕРКА ПУТИ: убедимся, что мы смотрим в правильное место
        path = f"data/{folder}/"
        if not os.path.exists(path):
            continue
            
        files = [f for f in os.listdir(path) if f.endswith('.maf.gz')]
        
        for file_name in files:
            file_path = os.path.join(path, file_name)
            try:
                with gzip.open(file_path, 'rt') as f:
                    # Пропускаем комментарии до заголовка
                    header = None
                    for line in f:
                        if line.startswith('Hugo_Symbol'):
                            header = line.strip().split('\t')
                            break
                    if not header: continue
                    
                    # Читаем данные
                    df = pd.read_csv(f, sep='\t', names=header, comment='#', low_memory=False)
                    
                    # Ищем наши гены
                    found = df[df['Hugo_Symbol'].isin(target_genes)]
                    if not found.empty:
                        for _, mut in found.iterrows():
                            print(f"Bingo! Found {mut['Hugo_Symbol']} mutation in {case}")
                            mutations.append({
                                'Case_ID': case,
                                'Cancer': cancer,
                                'Gene': mut['Hugo_Symbol'],
                                'Variant_Class': mut.get('Variant_Classification', 'Unknown'),
                                'HGVSp': mut.get('HGVSp_Short', 'N/A')
                            })
            except Exception as e:
                print(f"Error reading {file_name}: {e}")

    # СОЗДАЕМ ФАЙЛ В ЛЮБОМ СЛУЧАЕ
    mut_df = pd.DataFrame(mutations)
    os.makedirs('results', exist_ok=True)
    mut_df.to_csv('results/detailed_mutations.tsv', sep='\t', index=False)
    print(f"\nDone. Found {len(mutations)} mutations. File saved to results/detailed_mutations.tsv")

if __name__ == "__main__":
    main()
