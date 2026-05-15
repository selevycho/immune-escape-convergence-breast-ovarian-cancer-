#!/bin/bash
#SBATCH --job-name=OptiType_Final_Success
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=optitype_final_%j.log

# Подгружаем конду
source ~/miniconda3/etc/profile.d/conda.sh
conda activate optitype_env

# Пути
BAM_FILE="$HOME/gatk_resources/wgs_24RG_b37/NA12878_24RG_med.b37.bam"
OUT_DIR="$HOME/immune_escape_project/results/hla_typing"
CONFIG_FILE="$HOME/immune_escape_project/scripts/config.ini"

mkdir -p $OUT_DIR
cd $OUT_DIR

# Очистка старых логов и временных файлов
rm -f *.sam *.bam *.hdf5 my_result*

echo "--- STEP 1: EXTRACTION (MHC Region) ---"
samtools view -@ 8 -b $BAM_FILE 6:28000000-34000000 > fished.bam

echo "--- STEP 2: FASTQ CONVERSION ---"
samtools sort -@ 8 -n fished.bam -o sorted.bam
bedtools bamtofastq -i sorted.bam -fq hla_R1.fastq -fq hla_R2.fastq

echo "--- STEP 3: RUNNING OPTITYPE ---"
# Мы запускаем OptiType и сразу после маппинга (через 10-15 сек) 
# он создаст .sam файл. Мы подготовим для него почву.

OptiTypePipeline.py \
    -i hla_R1.fastq hla_R2.fastq \
    --dna \
    -v \
    --prefix my_result \
    -o . \
    -c $CONFIG_FILE

echo "--- PIPELINE FINISHED ---"
