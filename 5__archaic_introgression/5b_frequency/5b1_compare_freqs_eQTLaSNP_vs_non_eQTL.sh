#!/bin/bash
################################################################################
################################################################################
# File name: 5b1_compare_freqs_eQTLaSNP_vs_non_eQTL.sh
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compare frequency of (pruned) aSNPs between eQTL and non eQTL
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=5b1
#SBATCH --output=${LOG_DIR}/5b1__%J.log
#SBATCH -J compareFreq
#SBATCH --mem=8G
################################################################################
################################################################################
# Command history
# number of sets in /3__eQTL_mapping/SumStats/All_eQTL_snpsSets.txt.gz',
# NSET=198
# sbatch --array=1-${NSET} ./5b1_compare_freqs_eQTLaSNP_vs_non_eQTL.sh CEU
# sbatch --array=1-${NSET} ./5b1_compare_freqs_eQTLaSNP_vs_non_eQTL.sh CHS

# Setup
## libs
module purge
module load R/4.1.0
module load tabix
module load samtools/1.10

## args
num_test=$SLURM_ARRAY_TASK_ID
pop=$1

## wd
echo "Starting time: `date`"

echo ${pop} ${num_test}
## run R
  Rscript ./5b1_compare_freqs_eQTLaSNP_vs_non_eQTL.R ${num_test} ${pop}

echo "End time: `date`"
