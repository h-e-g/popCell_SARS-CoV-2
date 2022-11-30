#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mem 100G
#SBATCH -o "1d1_%J.log"
#SBATCH -J covariates

################################################################################
# Commands: on expression

# Compute covariates per lineage, without cell type proportions, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv 0 --prop FALSE --run_id lineage_condition__noCellProp_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv -1 --prop FALSE --run_id lineage_condition__noCellProp_SVs

# Compute covariates per lineage, with intra-lineage cell type proportions, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv 0 --prop TRUE --run_id lineage_condition__CellPropLineage_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv -1 --prop TRUE --run_id lineage_condition__CellPropLineage_SVs

# Compute covariates per cell type, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype celltype --sv 0 --prop FALSE --run_id celltype_condition__noCellProp_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype celltype --sv -1 --prop FALSE --run_id celltype_condition__noCellProp_SVs

################################################################################

################################################################################
# Commands: on responses

# Compute covariates per lineage, without cell type proportions, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv 0 --prop FALSE --logfc TRUE --run_id lineage_condition_logFC__noCellProp_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv -1 --prop FALSE --logfc TRUE --run_id lineage_condition_logFC__noCellProp_SVs

# Compute covariates per lineage, with intra-lineage cell type proportions, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv 0 --prop TRUE  --logfc TRUE --run_id lineage_condition_logFC__CellPropLineage_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype lineage --sv -1 --prop TRUE --logfc TRUE --run_id lineage_condition_logFC__CellPropLineage_noSVs

# Compute covariates per cell type, with or without SVs
# sbatch ./1d1_computeCovariates.sh --celltype celltype --sv 0 --prop FALSE --logfc TRUE --run_id celltype_condition_logFC__noCellProp_noSVs
# sbatch ./1d1_computeCovariates.sh --celltype celltype --sv -1 --prop FALSE --logfc TRUE --run_id celltype_condition_logFC__noCellProp_noSVs

################################################################################


# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="1__transcriptome_processing/1d__covariate_computation"

Rscript ${SCRIPT_DIR}/1d1_computeCovariates.R $*

exit
