#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 250G
#SBATCH -o "log/08_%J.log"
#SBATCH -J pDR_3
#SBATCH  --array=1-22

###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
VCF_FILE="file_cleaned_merged.vcf.gz"
OUT_DIR="${EIP}/project/pop_eQTL/data/4_natural_selection/clues/phased_snps"
POPCELL_DATA="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"

module load plink
module load R
module load samtools
module load tabix

## index the files
bcftools index ${VCF_FILE} -t

###############################################################################
                                 # Phasing #
###############################################################################

### Software:
  # Shapeit4: already a module on Maestro
  # Download map file (already done by Maxime): 
  # cd ${EIP}/resources/references
  # wget https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz
    # gunzip genetic_maps.b38.tar.gz
    # tar -xvf genetic_maps.b38.tar


#module load shapeit/r837
module load plink/1.90b5.4
module load tabix
module load samtools/1.10
module load vcftools/0.1.16
module load shapeit4/4.2.1

###############################################################################

# check with Maxime how he performed the phasing

srun -p geh  shapeit4 --input ${VCF_FILE}  --map ${EIP}/resources/references/genetic_maps.b38/chr${SLURM_ARRAY_TASK_ID}.b38.gmap.gz --region $SLURM_ARRAY_TASK_ID --output ${OUT_DIR}/Geno_b38_473Ind_3723840snps_chr${SLURM_ARRAY_TASK_ID}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --log ${OUT_DIR}/phased.log

###############################################################################

module load samtools

cd ${OUT_DIR}

srun -p geh bcftools index file_${SLURM_ARRAY_TASK_ID}_phased.vcf.gz -t

