#!/bin/bash
#SBATCH --job-name=bwa_mem
#SBATCH --partition=cpu-single       # На bwUniCluster 3.0 используем эту очередь
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # 16 ядер ускорят процесс в 16 раз
#SBATCH --mem=64G                    # Памяти с запасом для сортировки
#SBATCH --time=08:00:00              # Даем 8 часов, чтобы точно не вылетело
#SBATCH --output=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/bwa_%j.log
#SBATCH --error=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/bwa_%j.err

# Пути
BASE_DIR="/home/fr/fr_fr/fr_os136/immune_escape_project"
IN_DIR="${BASE_DIR}/$1"
OUT_DIR="${BASE_DIR}/$2"
SAMPLE=$3
REF="${BASE_DIR}/data/ref/hg38.fa"

# Активация среды
source /home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
conda activate bio_work

mkdir -p $OUT_DIR

echo "--- Start Alignment: $SAMPLE ---"

# 1. BWA MEM -> BAM (с Read Groups для Mutect2)
bwa mem -t 16 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" $REF \
    ${IN_DIR}/*_1.fastq.gz ${IN_DIR}/*_2.fastq.gz | \
    samtools view -@ 4 -Sb - > ${OUT_DIR}/${SAMPLE}_unsorted.bam

echo "--- Start Sorting: $SAMPLE ---"
# 2. Сортировка (необходима для поиска мутаций)
samtools sort -@ 16 ${OUT_DIR}/${SAMPLE}_unsorted.bam -o ${OUT_DIR}/${SAMPLE}_sorted.bam

echo "--- Start Indexing BAM: $SAMPLE ---"
# 3. Индексация BAM (создает .bai)
samtools index ${OUT_DIR}/${SAMPLE}_sorted.bam

# Удаляем временный файл
rm ${OUT_DIR}/${SAMPLE}_unsorted.bam

echo "--- Finished: $SAMPLE ---"
