#!/bin/bash
################################################################################
################################################################################
# File name: 0b1__demuxlet.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-library genotype-based barcode demultiplexing
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
LIB_COUNT=134

################################################################################
# SLURM options
#SBATCH --job-name=0b1
#SBATCH --output=${LOG_DIR}/0b1__%A_%a.log
#SBATCH --mem=128G
#SBATCH --ntasks=16
#SBATCH --array=1-${LIB_COUNT}
################################################################################

LIBID=${SLURM_ARRAY_TASK_ID}

################################################################################
# Setup

TASKDIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
PROJECT="HN00163853"

# load modules
module load popscle/0.1-beta

# create temporary VCFs, run once
#cd /pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data
#
#export PATH=/pasteur/sonic/homes/yaaquino/miniconda3/bin:$PATH
#
#source activate gatk
#
#mkdir liftover_tmp
#
#/pasteur/zeus/projets/p02/IGSR/Automated_Pipeline_For_WGS/Resources/Software/bin/gatk-4.1.2.0/gatk LiftoverVcf --I EvoImmunoPop_Omni5_473x3723480_chrom.vcf -O EvoImmunoPop_Omni5_473x3723480_lifted38.vcf --CHAIN filter_chain.txt --REJECT rejected_variants.vcf -R human_iav_sars/genome.fa --TMP_DIR liftover_tmp


# create library-specific VCF
#/pasteur/zeus/projets/p02/IGSR/Automated_Pipeline_For_WGS/Resources/Software/bin/bcftools-1.11/bcftools view --force-samples -S individuals_in_L${LIBID}_dmx.tsv EvoImmunoPop_Omni5_473x3723480_lifted38.vcf > EvoImmunoPop_Omni5_473x3723480.vcf


VCF="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/EvoImmunoPop_Omni5_473x3723480_lifted38.vcf"

################################################################################
# Run Demuxlet

mkdir -p ${TASKDIR}/demultiplex/${PROJECT}/demuxlet/L${LIBID}

echo "demuxlet start."

popscle demuxlet --sam ${TASKDIR}/aligned/${PROJECT}/L${LIBID}.Aligned.sortedByCoord.out.bam --vcf ${VCF} --field GT --group-list ${TASKDIR}/aligned/${PROJECT}/L${LIBID}.Solo.out/Gene/filtered/barcodes.tsv --out ${TASKDIR}/demultiplex/${PROJECT}/demuxlet/L${LIBID}/L${LIBID}_demuxlet

echo "demuxlet stop."
