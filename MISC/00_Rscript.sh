#!/bin/bash
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END

module purge # remove modules that may have been loaded by mistake
module load R/4.1.0
module load samtools # needed to call bcftools from within R (see ./MISC/querySNPs.R)
module load SLiM/4.0.1 # needed to call SLiM from within R (see ../MISC/misc_SLiM.R)

SCRIPT_DIR=""
SCRIPT=$1
shift;

Rscript ${SCRIPT_DIR}/${SCRIPT} $* ${SLURM_ARRAY_TASK_ID}
exit
