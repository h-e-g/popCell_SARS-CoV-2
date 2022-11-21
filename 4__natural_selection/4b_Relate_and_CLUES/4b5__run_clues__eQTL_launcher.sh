#!/bin/bash
################################################################################
################################################################################
# File name: 4b5__run_clues__eQTL_launcher.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Resample branch lengths per eQTL and run CLUES
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=4b5
#SBATCH --output=${LOG_DIR}/4b5__%A_%a.log
#SBATCH --mem=8G
#SBATCH --array=1-12753

################################################################################
# Setup

POP=$1
OUT_DIR="4__natural_selection/data/CLUES"

NUM=${SLURM_ARRAY_TASK_ID}

RSID=$(sed "${NUM}q;d" ${OUT_DIR}/eQTL_position.tsv | cut  -f1)
CHR=$(sed "${NUM}q;d" ${OUT_DIR}/eQTL_position.tsv | cut  -f2)
POS=$(sed "${NUM}q;d" ${OUT_DIR}/eQTL_position.tsv | cut  -f3)

CUTOFF=2000

################################################################################
# Command history

# sbatch ./4b5__run_clues__eQTL_launcher.sh CHS
# sbatch ./4b5__run_clues__eQTL_launcher.sh CEU
# sbatch ./4b5__run_clues__eQTL_launcher.sh YRI

################################################################################
# Run bash scripts

sh ./4b5__run_clues__eQTL.sh ${RSID} ${CHR} ${POS} ${POP} ${CUTOFF} eQTL
