#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#SBATCH -o /pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/logs/slurm-%A_%a_PBS.out
#SBATCH -J resamp
#SBATCH --mem=16G
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END

## libs
module purge
module load R/4.1.0
module load plink/1.90b6.16
module load tabix
module load samtools/1.10
module load vcftools/0.1.16
module load java/13.0.2
module load shapeit4/4.2.1

## args
num_test=$SLURM_ARRAY_TASK_ID
num_resamples=$1
pop=$2

## paths
SS="/pasteur/zeus/projets/p02/evoceania/Javier/Softwares_scripts/"
IGSR="/pasteur/zeus/projets/p02/IGSR/"
REFhg38="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/references/RNA/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
SP="/pasteur/appa/scratch/Javier_tmp/"

## wd
cd "/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier"

echo "Starting time: `date`"

## run R
Rscript ./scripts/0064_compute_enrichment_PBS_3factors.R ${num_test} ${num_resamples} ${pop}

echo "End time: `date`"
