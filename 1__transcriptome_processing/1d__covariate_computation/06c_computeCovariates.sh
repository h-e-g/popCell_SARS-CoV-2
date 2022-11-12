#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mem 100G
#SBATCH --mail-type=END
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH -o "log/06c_%J.log"
#SBATCH -J cov

################################################################################
# Commands: on expression

# Compute covariates per cell type, include intra-lineage cell type proportions or not
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop FALSE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop FALSE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop TRUE --proptype lineage
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop TRUE --proptype lineage

# Compute covariates per cell type, include intra- and extra-lineage cell type proportions
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop TRUE --proptype celltype
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop TRUE --proptype celltype

# Compute covariates per lineage, include intra-lineage cell type proportions or not
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop FALSE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop FALSE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop TRUE --proptype lineage
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop TRUE --proptype lineage

# Compute covariates per lineage, include intra- and extra-lineage cell type proportions
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop TRUE --proptype celltype
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop TRUE --proptype celltype
################################################################################

################################################################################
# Commands: on responses

# Compute covariates per cell type, include intra-lineage cell type proportions or not
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop FALSE --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop FALSE --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop TRUE --proptype lineage --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop TRUE --proptype lineage --logfc TRUE

# Compute covariates per cell type, include intra- and extra-lineage cell type proportions --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv 0 --subset FALSE --prop TRUE --proptype celltype --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype celltype --state condition --sv -1 --subset FALSE --prop TRUE --proptype celltype --logfc TRUE

# Compute covariates per lineage, include intra-lineage cell type proportions or not --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop FALSE --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop FALSE --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop TRUE --proptype lineage --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop TRUE --proptype lineage --logfc TRUE

# Compute covariates per lineage, include intra- and extra-lineage cell type proportions --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv 0 --subset FALSE --prop TRUE --proptype celltype --logfc TRUE
# sbatch ./06c_computeCovariates.sh --celltype lineage --state condition --sv -1 --subset FALSE --prop TRUE --proptype celltype --logfc TRUE
################################################################################ --logfc TRUE


# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/clean"

Rscript ${SCRIPT_DIR}/06c_computeCovariates.R $*

exit
