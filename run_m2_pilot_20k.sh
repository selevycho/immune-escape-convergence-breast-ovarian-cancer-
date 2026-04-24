#!/bin/bash
#SBATCH --job-name=Mutect2_Pilot_20k
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=mutect_pilot_%j.log

# ==============================================================================
# PROJECT: Immune Escape Convergence (Breast & Ovarian Cancer)
# STEP: Pilot Somatic Variant Calling Validation
# DATA: NA12878 (Downsampled WGS - 20k reads)
# REF: human_g1k_v37 (Chr 20, 21 only)
# 
# DESCRIPTION: This script validates the GATK4 environment and tests 
# interval-based calling to handle potential reference contig mismatches.
# ==============================================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

# Define paths
REF=$HOME/gatk_resources/mutect2_b37/human_g1k_v37.20.21.fasta
BAM=$HOME/gatk_resources/wgs_20k_b37/NA12878_20k.b37.bam

echo "Starting Mutect2 for chromosomes 20 and 21..."

# Step 1: Somatic Variant Calling
gatk Mutect2 \
  -R $REF \
  -I $BAM \
  -L 20 -L 21 \
  -O raw_variants.vcf

# Step 2: Variant Filtering
gatk FilterMutectCalls \
  -R $REF \
  -V raw_variants.vcf \
  -O filtered_variants.vcf

echo "Analysis complete. Output saved to filtered_variants.vcf"
echo "DONE!"
