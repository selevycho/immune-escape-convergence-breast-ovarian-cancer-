#!/bin/bash
#SBATCH --job-name=LOHHLA_RUN
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/logs/lohhla_%j.log

# Проверка: передал ли пользователь ID пациента (например, Patient2)
if [ -z "$1" ]; then
    echo "Ошибка: Нужно указать ID пациента. Пример: sbatch run_lohhla_generic.sh Patient2"
    exit 1
fi

PATIENT_ID=$1

# 1. Активация окружения
source /home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
conda activate lohhla_env

# 2. Определение базовых путей
BASE_DIR="/home/fr/fr_fr/fr_os136/immune_escape_project"
LOHHLA_DIR="${BASE_DIR}/soft/lohhla"
OUT_DIR="${BASE_DIR}/results/LOHHLA/${PATIENT_ID}"

mkdir -p $OUT_DIR
mkdir -p ${BASE_DIR}/scripts/logs

# Настройка путей к файлам (исходя из твоей структуры папок)
# ВАЖНО: Убедись, что для Patient3 BAM-файлы лежат в папках PEO... согласно логике
HLA_PATH="${BASE_DIR}/results/HLA/PEO1/${PATIENT_ID}_alleles.txt"
NORMAL_BAM="${BASE_DIR}/results/BAM/PEO1/PEO1_sorted.bam"
TUMOR_BAM="${BASE_DIR}/results/BAM/PEO4/PEO4_sorted.bam"

# 3. Запуск LOHHLA
Rscript ${LOHHLA_DIR}/LOHHLAscript.R \
  --patientId ${PATIENT_ID} \
  --outputDir ${OUT_DIR} \
  --normalBAMfile ${NORMAL_BAM} \
  --tumorBAMfile ${TUMOR_BAM} \
  --hlaPath ${HLA_PATH} \
  --HLAfastaLoc ${LOHHLA_DIR}/data/abc_complete.fasta \
  --HLAexonLoc ${LOHHLA_DIR}/data/hla.dat \
  --CopyNumLoc FALSE \
  --mappingStep TRUE \
  --minCoverageFilter 10 \
  --numMisMatch 1 \
  --plottingStep TRUE \
  --ignoreWarnings TRUE

echo "LOHHLA finished for ${PATIENT_ID}"
