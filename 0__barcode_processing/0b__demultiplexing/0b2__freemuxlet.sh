#!/bin/bash
################################################################################
################################################################################
# File name: 0b2__freemuxlet.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-library genotype-free barcode demultiplexing
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
LIB_COUNT=134

################################################################################
# SLURM options
#SBATCH --job-name=0b2
#SBATCH --output=${LOG_DIR}/0b2__%A_%a.log
#SBATCH --mem=128G
#SBATCH --ntasks=16
#SBATCH --array=1-${LIB_COUNT}
################################################################################

LIBID=${SLURM_ARRAY_TASK_ID}
VCF="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/EvoImmunoPop_Omni5_473x3723480_lifted38.vcf"

################################################################################
# Setup

TASKDIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
PROJECT="HN00163853"

# load modules
module load popscle/0.1-beta

nsample=$(wc -l individuals_in_L${LIBID}_dmx.tsv | cut -d" " -f1)

################################################################################
# Run Freemuxlet

mkdir -p ${TASKDIR}/demultiplex/${PROJECT}/freemuxlet/L${LIBID}

echo "dsc-pileup start."

popscle dsc-pileup --sam ${TASKDIR}/aligned/${PROJECT}/L${LIBID}.Aligned.sortedByCoord.out.bam --vcf ${VCF} --group-list ${TASKDIR}/aligned/${PROJECT}/L${LIBID}.Solo.out/Gene/filtered/barcodes.tsv --out ${TASKDIR}/demultiplex/${PROJECT}/freemuxlet/L${LIBID}/L${LIBID}_pileup

echo "dsc-pileup stop."

echo "freemuxlet start."

popscle freemuxlet --plp ${TASKDIR}/demultiplex/${PROJECT}/freemuxlet/L${LIBID}/L${LIBID}_pileup --out ${TASKDIR}/demultiplex/${PROJECT}/freemuxlet/L${LIBID}/L${LIBID}_freemuxlet --nsample ${nsample}

echo "freemuxlet stop."
