#!/bin/bash
################################################################################
################################################################################
# File name: 1c2__pseudobulk_batch_correction.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Pseudobulk computation from filtered SingleCellExperiment object
# Batch launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --mem 1000G
#SBATCH -o "${LOG_DIR}/1c2__%J.log"
#SBATCH -J 1c2
################################################################################

################################################################################
# Usage
# sbatch ./1c2__pseudobulk_batch correction.sh --celltype celltype --state condition
# sbatch ./1c2__pseudobulk_batch correction.sh --celltype lineage --state condition
################################################################################

################################################################################
# Setup

# load modules
module load R/4.1.0

SCRIPT_DIR="1__transcriptome_processing/1c__pseudobulk_computation"

Rscript ${SCRIPT_DIR}/1c2__pseudobulk_batch correction.R $*

exit
