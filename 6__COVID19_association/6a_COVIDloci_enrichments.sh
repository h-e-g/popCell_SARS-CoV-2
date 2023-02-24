#!/bin/bash
################################################################################
################################################################################
# File name: 6a_COVIDloci_enrichments.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: test enrichment of eQTL and reQTLS in COVID risk loci
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=6a
#SBATCH --output=${LOG_DIR}/6a__%J.log
#SBATCH -J CovidEnrich
#SBATCH --mem=120G
################################################################################
################################################################################
# Command history
# number of sets in /3__eQTL_mapping/SumStats/All_eQTL_snpsSets.txt.gz',
# NSET=198
# sbatch --array=1-${NSET} ./6a_COVIDloci_enrichments.sh

# Setup
## libs
module purge
module load R/4.1.0
module load tabix
module load samtools/1.10

## args
num_set=$SLURM_ARRAY_TASK_ID

## wd
echo "Starting time: `date`"

echo ${num_set}
## run R
  Rscript ./6a_COVIDloci_enrichments.R --nsamp 10000 --set ${num_set}

echo "End time: `date`"
