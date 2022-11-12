#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mem 1000G
#SBATCH --mail-type=END
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH -o "log/06b_%J.log"
#SBATCH -J bat_c
#SBATCH --dependency=afterok:5486237

################################################################################
# Command:
# sbatch ./06b_batch_correct.sh --celltype celltype --state condition
# sbatch ./06b_batch_correct.sh --celltype lineage --state condition
################################################################################

# load modules you need:
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/clean"

Rscript ${SCRIPT_DIR}/06b_batch_correct.R $*

exit
