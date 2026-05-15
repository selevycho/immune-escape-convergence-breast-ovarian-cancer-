import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import os

def generate_nice_dashboard(vcf_path):
    if not os.path.exists(vcf_path):
        vcf_path = os.path.join('..', vcf_path)

    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    df = df[df['FILTER'] == 'PASS'].copy()

    sample_col = df.columns[-1]
    df['VAF'] = df[sample_col].apply(lambda x: float(x.split(':')[2]) if len(x.split(':')) > 2 else 0)
    
    # Группируем замены для чистоты графика
    def simplify_sub(row):
        sub = f"{row['REF']}>{row['ALT']}"
        standard = ['C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G']
        # Инвертируем G>A в C>T для стандартного представления (онкологический стандарт)
        mapping = {'G>A':'C>T', 'G>C':'C>G', 'G>T':'C>A', 'A>T':'T>A', 'A>G':'T>C', 'A>C':'T>G'}
        return mapping.get(sub, sub if len(sub)==3 else 'Indels/Other')

    df['Standard_Sub'] = df.apply(simplify_sub, axis=1)

    # Рисуем
    plt.style.use('seaborn-v0_8-muted')
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # 1. VAF Plot
    sns.histplot(ax=axes[0], data=df, x='VAF', bins=30, kde=True, color='skyblue', label='Observed VAF')
    axes[0].axvline(0.5, color='red', linestyle='--', label='Germline Expectation (0.5)')
    axes[0].set_title('Variant Allele Frequency (VAF) Distribution', fontsize=14)
    axes[0].set_xlabel('VAF (Allele Fraction)', fontsize=12)
    axes[0].legend()

    # 2. Mutational Spectrum (Чистый)
    order = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'Indels/Other']
    sns.countplot(ax=axes[1], data=df, x='Standard_Sub', order=order, palette='viridis')
    axes[1].set_title('Mutational Spectrum (Standardized)', fontsize=14)
    axes[1].set_xlabel('Substitution Type', fontsize=12)
    
    plt.suptitle(f"Pilot Analysis Report: {len(df)} PASS Variants Found", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig('../nice_dashboard.png', dpi=300)
    print("�� Новый красивый дашборд готов: nice_dashboard.png")

generate_nice_dashboard('filtered_variants.vcf')
