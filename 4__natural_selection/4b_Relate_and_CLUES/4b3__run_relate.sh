#!/bin/bash
################################################################################
################################################################################
# File name: 4b3__run_relate.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Run Relate on filtered files
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b3
#SBATCH --output=${LOG_DIR}/4b2__%A_%a.log
#SBATCH --mem=50G
#SBATCH  --array=1-22

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

################################################################################
# Command history

# sbatch ./4b3__run_relate.sh CHS
# sbatch ./4b3__run_relate.sh CEU
# sbatch ./4b3__run_relate.sh YRI

################################################################################
# Run Relate

mkdir ${OUT_DIR}/${POP}
cd ${OUT_DIR}/${POP}

${PATH_RELATE}/bin/Relate --mode All \
  -m 1.25e-8 \
  -N 30000 \
  --haps ${IN_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}_input.haps.gz \
  --sample ${IN_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}_input.sample.gz \
  --map ${GMAP_DIR}/chr${CHR}.b38.gmap.gz \
  --seed 1 \
  -o relate_output_chr${CHR}
