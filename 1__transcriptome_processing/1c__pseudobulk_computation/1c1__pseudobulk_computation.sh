#!/bin/bash
################################################################################
################################################################################
# File name: 1c1__pseudobulk_computation.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Pseudobulk computation from filtered SingleCellExperiment object
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
LIB_COUNT=135

################################################################################
# SLURM options
#SBATCH --mem 250G
#SBATCH -o "${LOG_DIR}/1c1__%J.log"
#SBATCH -J 1c1
################################################################################

################################################################################
# Usage
# sbatch ./06a_compute_pseudoBulk_from_sce.sh --celltype celltype --state condition
# sbatch ./06a_compute_pseudoBulk_from_sce.sh --celltype lineage --state condition
################################################################################

################################################################################
# Setup

# load modules
module load R/4.1.0

SCRIPT_DIR="1__transcriptome_processing/1c__pseudobulk_computation"

Rscript ${SCRIPT_DIR}/1c1__pseudobulk_computation.R $*

exit
