#!/bin/bash
#SBATCH --job-name=hg38_idx
#SBATCH --partition=cpu-single            # <--- КЛЮЧЕВОЕ ИЗМЕНЕНИЕ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb                   # Увеличил до 64GB для стабильности на bwUniCluster
#SBATCH --time=02:00:00              # Теперь Slurm не выкинет тебя через 30 минут
#SBATCH --output=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/bwa_idx_%j.log
#SBATCH --error=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/bwa_idx_%j.err

# Пути
REF_DIR="/home/fr/fr_fr/fr_os136/immune_escape_project/data/ref"
CONDA_PATH="/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh"

source $CONDA_PATH
conda activate bio_work

cd $REF_DIR

echo "Cleaning old bits..."
rm -f hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa

echo "Starting BWA indexing..."
bwa index hg38.fa
echo "Indexing finished!"
