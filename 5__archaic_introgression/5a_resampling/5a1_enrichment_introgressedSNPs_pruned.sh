#!/bin/bash
################################################################################
################################################################################
# File name: 5a1_enrichment_introgressedSNPs_pruned.sh
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute excess of introgressed SNPs among various eQTL lists
# relative to random matched SNPs, and provide a p-value
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=5a1
#SBATCH --output=${LOG_DIR}/5a1__%J.log
#SBATCH -J resamp
#SBATCH --mem=40G
################################################################################
################################################################################
# Command history
# number of sets in /3__eQTL_mapping/SumStats/All_eQTL_snpsSets.txt.gz',
# NSET=198
# sbatch --array=1-${NSET} ./5a1_enrichment_introgressedSNPs_pruned.sh 10000 CEU
# sbatch --array=1-${NSET} ./5a1_enrichment_introgressedSNPs_pruned.sh 10000 CHS

# Setup
## libs
module purge
module load R/4.1.0
module load samtools/1.10

## args
num_test=$SLURM_ARRAY_TASK_ID  # number of the tested set
num_resamples=$1 # 10000
pop=$2

echo "Starting time: `date`"

## run R
  Rscript ./5a1_enrichment_introgressedSNPs_pruned.R ${num_test} ${num_resamples} ${pop}

echo "End time: `date`"
