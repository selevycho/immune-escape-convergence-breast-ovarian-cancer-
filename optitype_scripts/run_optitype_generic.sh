#!/bin/bash
#SBATCH --job-name=opti_gen
#SBATCH --partition=devel
#SBATCH --output=/home/fr/fr_fr/fr_os136/immune_escape_project/optitype_scripts/opti_%j.log
#SBATCH --error=/home/fr/fr_fr/fr_os136/immune_escape_project/optitype_scripts/opti_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00

# --- Проверка аргументов ---
if [ "$#" -ne 3 ]; then
    echo "Использование: sbatch run_optitype_generic.sh <ПУТЬ_ОТ_КОРНЯ_ПРОЕКТА> <ПАПКА_РЕЗУЛЬТАТОВ> <ИМЯ_ОБРАЗЦА>"
    exit 1
fi

# Используем полные пути, чтобы SLURM не терялся
BASE_DIR="/home/fr/fr_fr/fr_os136/immune_escape_project"
IN_DIR="${BASE_DIR}/$1"
OUT_BASE="${BASE_DIR}/$2"
PREFIX=$3

# Находим файлы
R1=$(ls ${IN_DIR}/*_1.fastq.gz)
R2=$(ls ${IN_DIR}/*_2.fastq.gz)

# Настройка среды
source /home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
conda activate optitype_stable
OPTITYPE_BIN="/home/fr/fr_fr/fr_os136/miniconda3/envs/optitype_stable/bin/OptiTypePipeline.py"

mkdir -p ${OUT_BASE}

# Запуск
python $OPTITYPE_BIN -i $R1 $R2 --dna -v -o ${OUT_BASE} --prefix ${PREFIX}
