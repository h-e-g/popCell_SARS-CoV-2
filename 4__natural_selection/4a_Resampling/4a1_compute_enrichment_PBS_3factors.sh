#!/bin/bash
################################################################################
################################################################################
# File name: 4a1_compute_enrichment_PBS_3factors.sh
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute excess of high PBS snps among various eQTL lists
# relative to random matched SNPs, and provide a p-value
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=4a1
#SBATCH --output=${LOG_DIR}/4a1__%J.log
#SBATCH -J resamp
#SBATCH --mem=16G
################################################################################
################################################################################
# Command history
# number of sets in /3__eQTL_mapping/SumStats/All_eQTL_snpsSets.txt.gz',
# NSET=198
# sbatch --array=1-${NSET} ./4a1_compute_enrichment_PBS_3factors.sh 10000 CEU
# sbatch --array=1-${NSET} ./4a1_compute_enrichment_PBS_3factors.sh 10000 CHS
# sbatch --array=1-${NSET} ./4a1_compute_enrichment_PBS_3factors.sh 10000 YRI

# Setup

# load modules
module purge
module load R/4.1.0
module load samtools/1.10

## args
num_test=$SLURM_ARRAY_TASK_ID
num_resamples=$1
pop=$2

echo "Starting time: `date`"

## run R
Rscript ./4a1_compute_enrichment_PBS_3factors.R ${num_test} ${num_resamples} ${pop}

echo "End time: `date`"
