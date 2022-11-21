#!/bin/bash
################################################################################
################################################################################
# File name: 4b1__convert_to_relate.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Convert One Thousand Genome Project's v3 VCFs to Relate input files
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b1
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

POP=$1
CHR=${SLURM_ARRAY_TASK_ID}
VCF_DIR="4__natural_selection/data/1KG/1KG_VCF"
VCF_PREFIX="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants_BiallelicSNPs_PASS_CR95_HWE_MER"
VCF_IN="${VCF_PREFIX}_shapeit4seq_8thread.vcf.gz"
VCF_OUT="${VCF_PREFIX}_shapeit4seq_8thread_${POP}.vcf.gz"
IN_DIR="4__natural_selection/data/1KG/1KG_VCF/phased"
OUT_DIR="4__natural_selection/data/1KG/1KG_VCF/RELATE"
META_DIR="4__natural_selection/data/1KG/1KG_Metadata"

# keep only unrelated individuals in each population
cd ${IN_DIR}
mkdir ${IN_DIR}/${POP}

# extract IDs
awk '{if (NR==FNR) {POP_SAMPLES[$1]=$1} else { if ($1 in POP_SAMPLES) {print $0}}}' \
  <(awk -F"\t" -v MYPOP=${POP} '{if ($4==MYPOP) print $1}' ${META_DIR}/igsr_samples.tsv) \
  ${META_DIR}/1kGP.2504_independent_samples.pedigree_info.txt > ${POP}/1KG_${POP}_IND.tsv

# subset phased VCFs
cd ${IN_DIR}/${POP}

${PATH_BCFTOOLS}/bcftools view -S 1KG_${POP}_IND.tsv -o ${VCF_OUT} ${IN_DIR}/${VCF_IN}

${PATH_BCFTOOLS}/bcftools index ${VCF_OUT} -t

################################################################################
# Command history

# sbatch ./4b1__convert_to_relate.sh CHS
# sbatch ./4b1__convert_to_relate.sh CEU
# sbatch ./4b1__convert_to_relate.sh YRI

################################################################################
# Convert VCFs

cd ${OUT_DIR}

${PATH_RELATE}/RelateFileFormats --mode ConvertFromVcf \
  --haps "${OUT_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}.haps" \
  --sample "${OUT_DIR}/${VCF_PREFIX}_shapeit4seq_8thread_${POP}.sample" \
  -i ${IN_DIR}/${POP}/${VCF_OUT}
