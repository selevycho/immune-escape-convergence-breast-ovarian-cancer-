#!/bin/bash
#SBATCH --job-name=mutect2
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/mutect2_%j.log

# Пути (проверь, чтобы они соответствовали твоей структуре)
BASE_DIR="/home/fr/fr_fr/fr_os136/immune_escape_project"
REF="${BASE_DIR}/data/ref/hg38.fa"
BAM_DIR="${BASE_DIR}/results/BAM"
OUT_DIR="${BASE_DIR}/results/VCF"

TUMOR_BAM=$1
NORMAL_BAM=$2
TUMOR_NAME=$3
NORMAL_NAME=$4
OUT_PREFIX=$5

# Активация среды
source /home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
conda activate bio_work

mkdir -p $OUT_DIR

echo "Running Mutect2: Tumor($TUMOR_NAME) vs Normal($NORMAL_NAME)"

gatk Mutect2 \
    -R $REF \
    -I ${BAM_DIR}/${TUMOR_NAME}/${TUMOR_BAM} \
    -I ${BAM_DIR}/${NORMAL_NAME}/${NORMAL_BAM} \
    -normal $NORMAL_NAME \
    -O ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
