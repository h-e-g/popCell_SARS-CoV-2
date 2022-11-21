#!/bin/bash
################################################################################
################################################################################
# File name: 4b6__run_clues__reQTL_launcher.sh
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
#SBATCH --output=${LOG_DIR}/4b6__%A_%a.log
#SBATCH --mem=8G
#SBATCH --array=1-1505

################################################################################
# Setup

POP=$1
OUT_DIR="4__natural_selection/data/CLUES"

NUM=${SLURM_ARRAY_TASK_ID}

RSID=$(sed "${NUM}q;d" ${OUT_DIR}/reQTL_position.tsv | cut  -f1)
CHR=$(sed "${NUM}q;d" ${OUT_DIR}/reQTL_position.tsv | cut  -f2)
POS=$(sed "${NUM}q;d" ${OUT_DIR}/reQTL_position.tsv | cut  -f3)

CUTOFF=2000

################################################################################
# Command history

# sbatch ./4b5__run_clues__reQTL_launcher.sh CHS
# sbatch ./4b5__run_clues__reQTL_launcher.sh CEU
# sbatch ./4b5__run_clues__reQTL_launcher.sh YRI

################################################################################
# Run bash scripts

sh ./4b6__run_clues__reQTL.sh ${RSID} ${CHR} ${POS} ${POP} ${CUTOFF} eQTL
