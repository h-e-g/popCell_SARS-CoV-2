#!/bin/bash
################################################################################
################################################################################
# File name: 1a1__quality_control.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-library and library-aggregated transcriptome quality control
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
LIB_COUNT=134

################################################################################
# SLURM options
#SBATCH --job-name=1a1
#SBATCH --output=${LOG_DIR}/1a1__%A_%a.log
#SBATCH --mem=100G
#SBATCH --array=1-${LIB_COUNT}
################################################################################

LIB=${SLURM_ARRAY_TASK_ID}

################################################################################
# Setup

# load modules
module load R/4.1.0

################################################################################
# Run R script

echo "Quality control for barcodes from library ${LIB} of ${LIB_COUNT}"

Rscript ./1a1__quality_control__per_library.R ${LIB}
