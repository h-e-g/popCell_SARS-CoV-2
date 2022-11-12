#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 50G
#SBATCH -o "log/%J_prepare_relate.log"
#SBATCH -J prp
#SBATCH  --array=1-22

# sbatch ./16c_prepare_relate.sh CHS
# sbatch ./16c_prepare_relate.sh CEU
# sbatch ./16c_prepare_relate.sh YRI

###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38"
OUT_DIR="${EIP}/project/pop_eQTL/data/4_natural_selection/clues/phased_snps"
POPCELL_DATA="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"

module load plink
module load R
module load samtools
module load tabix
module load vcftools/0.1.16

################################################################################

# prepare input files
PATH_RELATE="/pasteur/zeus/projets/p02/IGSR/Software/relate_v1.1.8_x86_64_static"
VCF_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05_chr"
PATH_ANCESTOR="/pasteur/zeus/projets/p02/IGSR/Human_Ancestor/GRCh38/homo_sapiens_ancestor_GRCh38"
PATH_MASK="/pasteur/zeus/projets/p02/IGSR/Genome_Masks/1KG/b38/PilotMask"

POP=$1

mkdir ${INPUT_DIR}/${POP}/relate_input
cd ${INPUT_DIR}/${POP}/relate_input

chr=${SLURM_ARRAY_TASK_ID}

${PATH_RELATE}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
                --haps ${INPUT_DIR}/${POP}/${VCF_PREFIX}${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz.haps \
                --sample ${INPUT_DIR}/${POP}/${VCF_PREFIX}${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz.sample \
                --ancestor ${PATH_ANCESTOR}/homo_sapiens_ancestor_${chr}.fa \
                --mask ${PATH_MASK}/20160622.chr${chr}.pilot_mask.fasta.gz \
                -o ${VCF_PREFIX}${chr}.filtered.shapeit2-duohmm-phased_${POP}_input
