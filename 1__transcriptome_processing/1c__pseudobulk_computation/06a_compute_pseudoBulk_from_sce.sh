#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mem 250G
#SBATCH --mail-type=END
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH -o "log/06a_%J.log"
#SBATCH -J pb_l

################################################################################
# Command:
# sbatch ./06a_compute_pseudoBulk_from_sce.sh --celltype celltype --state condition 
# sbatch ./06a_compute_pseudoBulk_from_sce.sh --celltype lineage --state condition 
################################################################################

# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/clean"

Rscript ${SCRIPT_DIR}/06a_compute_pseudoBulk_from_sce.R $*

exit
