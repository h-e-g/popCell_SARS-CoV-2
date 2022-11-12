#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p geh
#SBATCH --mem 250G
#SBATCH -o "log/08_AFBEUB_%A.log"
#SBATCH -J pDR_AE3
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END

################################################################################
# popDEGs

# Not adjust on cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_ --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_ --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_ --cellprop FALSE --perm 1 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_ --cellprop FALSE --perm 1 --runid 220409

# Adjust on intra-lineage cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage --cellprop FALSE --perm 1 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage --cellprop FALSE --perm 1 --runid 220409

# Adjust on intra- and extra-lineage cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition__CellPropCelltype --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition__CellPropCelltype --cellprop FALSE --perm 0 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition__CellPropCelltype --cellprop FALSE --perm 1 --runid 220409
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition__CellPropCelltype --cellprop FALSE --perm 1 --runid 220409

################################################################################
# popDRGs

# Not adjust on cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC_ --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC_ --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC_ --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC_ --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE

# Adjust on intra-lineage cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC__CellPropLineage --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC__CellPropLineage --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC__CellPropLineage --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC__CellPropLineage --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE

# Adjust on intra- and extra-lineage cell composition
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC__CellPropCelltype --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC__CellPropCelltype --cellprop FALSE --perm 0 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype lineage --state condition --covname lineage_condition_logFC__CellPropCelltype --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE
# sbatch ./08_PopulationDifference__AFB_EUB.sh --celltype celltype --state condition --covname celltype_condition_logFC__CellPropCelltype --cellprop FALSE --perm 1 --runid 220409 --logfc TRUE

# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/clean"

Rscript ${SCRIPT_DIR}/08_PopulationDifference__AFB_EUB.R $*

exit
