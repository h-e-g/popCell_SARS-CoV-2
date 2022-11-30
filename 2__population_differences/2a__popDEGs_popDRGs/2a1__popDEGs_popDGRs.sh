#!/bin/bash
################################################################################
################################################################################
# File name: 2a1__popDEGs_popDRGs.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use linear models to estimate population effects on immune response
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=2a1
#SBATCH --output=${LOG_DIR}/2a1__%J.log
#SBATCH --mem=250G
################################################################################

################################################################################
# Setup

# load modules
module load R/4.1.0

################################################################################
# Command history

# all arguments are described in the effector R script

# raw population expression differences

# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition__noCellProps_noSVs --perm 0 --runid 220409
# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition__noCellProps_noSVs --perm 1 --runid 220409

# cellular composition-adjusted population expression differences

# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition__CellPropLineage_noSVs --perm 0 --runid 220409
# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition__CellPropLineage_noSVs --perm 1 --runid 220409

# raw population response differences

# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition_logFC__noCellProps_noSVs --perm 0 --runid 220409 --logfc TRUE
# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition_logFC__noCellProps_noSVs --perm 1 --runid 220409 --logfc TRUE

# cellular composition-adjusted population response differences

# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition_logFC__CellPropLineage_noSVs --perm 0 --runid 220409 --logfc TRUE
# sbatch ./2a1__popDEGs_popDRGs.sh --celltype lineage --covname lineage_condition_logFC__CellPropLineage_noSVs --perm 1 --runid 220409 --logfc TRUE

################################################################################
# Run R script

Rscript ./2a1__popDEGs_popDRGs.R $*

exit
