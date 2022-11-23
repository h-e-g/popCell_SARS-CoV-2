#!/bin/bash
################################################################################
################################################################################
# File name: 4b6__run_clues__reQTL.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Resample branch lengths per reQTL and run CLUES
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b6
#SBATCH --output=${LOG_DIR}/4b6__%J.log
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

GMAP_DIR="4__natural_selection/data/1KG/1KG_GEN_MAP"
VCF_PREFIX="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants_BiallelicSNPs_PASS_CR95_HWE_MER"
IN_DIR="4__natural_selection/data/RELATE"
OUT_DIR="4__natural_selection/data/CLUES"
META_DIR="4__natural_selection/data/1KG/1KG_Metadata"

RSID=$1
CHR=$2
POS=$3
POP=$4
CUTOFF=$5

################################################################################
# Resample branch lengths

echo "Re-sampling branch lengths. ${CHR}:${POS} (${RSID}) in ${POP} (reQTL, CutOff=${CUTOFF})."

if test -f "${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/reQTL/${RSID}/relate_samplebranch_${RSID}.timeb"; then
    echo "relate_samplebranch_${RSID}.timeb  exists."
else
  cd ${IN_DIR}/${POP}
  mkdir -p ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/reQTL/${RSID}

  ${PATH_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
          -i relate_output_${CHR} \
          -o ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/eQTL/${RSID}/relate_samplebranch_${RSID} \
          -m 1.25e-8 \
          --coal relate_popsize_${POP}.coal
          --format b \
          --first_bp ${POS} \
          --last_bp ${POS} \
          --num_samples 100
fi

cd ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/reQTL/${RSID}

################################################################################
# Run CLUES

echo "Running CLUES. ${CHR}:${POS} (${RSID})."

source /pasteur/zeus/projets/p02/IGSR/Software/clues_conda/bin/activate

python3 ${PATH_CLUES}/inference_modified_path.py --times relate_output_${RSID} \
  --coal ${IN_DIR}/${POP}/relate_popsize_${POP}.coal \
  --out clues_output_${RSID} \
  --tCutoff ${CUTOFF}
