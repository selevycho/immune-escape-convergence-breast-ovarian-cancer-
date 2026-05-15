import pandas as pd
import io
import os

def inspect_vcf(vcf_path):
    # Проверка: если файла нет по указанному пути, попробуем на уровень выше
    if not os.path.exists(vcf_path):
        alt_path = os.path.join('..', vcf_path)
        if os.path.exists(alt_path):
            vcf_path = alt_path
        else:
            print(f"❌ Ошибка: Файл {vcf_path} не найден ни тут, ни в папке выше.")
            return

    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    pass_variants = df[df['FILTER'] == 'PASS']
    
    print("="*40)
    print(f"📊 ОТЧЕТ ПО ФАЙЛУ: {vcf_path}")
    print("="*40)
    print(f"✅ Всего надежных мутаций (PASS): {len(pass_variants)}")
    
    if len(pass_variants) > 0:
        subs = pass_variants['REF'] + ">" + pass_variants['ALT']
        print("\n🧬 Топ типов замен:")
        print(subs.value_counts().head(5))
    
    print("="*40)

inspect_vcf('filtered_variants.vcf')
