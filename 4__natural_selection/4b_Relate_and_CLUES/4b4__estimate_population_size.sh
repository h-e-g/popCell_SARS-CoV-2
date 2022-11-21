#!/bin/bash
################################################################################
################################################################################
# File name: 4b4__estimate_population_size.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Estimate population sizes
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b4
#SBATCH --output=${LOG_DIR}/4b4__%J.log
#SBATCH --mem=50G

################################################################################
# Setup

# load modules
module load plink
module load R
module load samtools
module load tabix
module load vcftools/0.1.16

PATH_BCFTOOLS="/pasteur/zeus/projets/p02/IGSR/Software/bcftools-1.15"
PATH_RELATE="/pasteur/zeus/projets/p02/IGSR/Software/relate_v1.1.8_x86_64_static/bin"
PATH_ANCESTOR="/pasteur/zeus/projets/p02/IGSR/Human_Ancestor/GRCh38/homo_sapiens_ancestor_GRCh38"
PATH_MASK="/pasteur/zeus/projets/p02/IGSR/Genome_Masks/1KG/b38/PilotMask"

POP=$1
CHR=${SLURM_ARRAY_TASK_ID}
GMAP_DIR="4__natural_selection/data/1KG/1KG_GEN_MAP"
VCF_PREFIX="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants_BiallelicSNPs_PASS_CR95_HWE_MER"
IN_DIR="4__natural_selection/data/1KG/1KG_VCF/RELATE/prepared"
OUT_DIR="4__natural_selection/data/RELATE"
META_DIR="4__natural_selection/data/1KG/1KG_Metadata"

# output poplabels file
cd ${IN_DIR}

awk -F"\t" '{if (NR==FNR) {POP_SAMPLES[$1]=$1} else { if ($1 in POP_SAMPLES) {print $1" "$4" "$6}}}' \
  "4__natural_selection/data/1KG/1KG_VCF/phased/${POP}/1KG_${POP}_IND.tsv" \
  ${META_DIR}/igsr_samples.tsv > poplabels_${POP}.tsv

################################################################################
# Command history

# sbatch ./4b4__estimate_population_size.sh CHS
# sbatch ./4b4__estimate_population_size.sh CEU
# sbatch ./4b4__estimate_population_size.sh YRI

################################################################################
# Estimate population sizes

cd ${OUT_DIR}/${POP}

${PATH_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i relate_output \
              -m 1.25e-8 \
              --seed 1 \
              --poplabels ${IN_DIR}/poplabels_${POP}.txt \
              -o relate_popsize_${POP} \
              --first_chr 1 \
              --last_chr 22
