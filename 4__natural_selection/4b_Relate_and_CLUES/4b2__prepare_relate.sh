#!/bin/bash
################################################################################
################################################################################
# File name: 4b2__prepare_relate.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Filter Relate input files
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b2
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
VCF_PREFIX="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants_BiallelicSNPs_PASS_CR95_HWE_MER"
IN_DIR="4__natural_selection/data/1KG/1KG_VCF/RELATE"
OUT_DIR="4__natural_selection/data/1KG/1KG_VCF/RELATE/prepared"
META_DIR="4__natural_selection/data/1KG/1KG_Metadata"

################################################################################
# Command history

# sbatch ./4b2__prepare_relate.sh CHS
# sbatch ./4b2__prepare_relate.sh CEU
# sbatch ./4b2__prepare_relate.sh YRI

################################################################################
# Prepare Relate input

mkdir ${OUT_DIR}
cd ${OUT_DIR}

${PATH_RELATE}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
  --haps ${IN_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}.haps \
  --sample ${IN_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}.sample \
  --ancestor ${PATH_ANCESTOR}/homo_sapiens_ancestor_${CHR}.fa \
  --mask ${PATH_MASK}/20160622.chr${CHR}.pilot_mask.fasta.gz \
  -o ${VCF_PREFIX}_shapeit4seq_8thread_${POP}_input
