#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mem 250G
#SBATCH --mail-type=END
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH -o "log/06d_%J.log"
#SBATCH -J pb_l

################################################################################
# Command:
# sbatch ./06d_compute_pseudoBulk_from_sce_nocelltype.sh --state condition
################################################################################

# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/clean"

Rscript ${SCRIPT_DIR}/06d_compute_pseudoBulk_from_sce_nocelltype.R $*

exit
