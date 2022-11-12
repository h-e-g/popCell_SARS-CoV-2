#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 50G
#SBATCH -o "log/%J_convert_to_relate.log"
#SBATCH -J cvt
#SBATCH  --array=1-22

# sbatch ./16b_convert_to_relate.sh CHS
# sbatch ./16b_convert_to_relate.sh CEU
# sbatch ./16b_convert_to_relate.sh YRI

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
module load vcftools/0.1.16

################################################################################

INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38"
META_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Metadata"

POP=$1

#cd ${INPUT_DIR}
#mkdir ${INPUT_DIR}/${POP}
#
## extract IDs of unrelated CHS individuals
#awk '{if (NR==FNR) {POP_SAMPLES[$1]=$1} else { if ($1 in POP_SAMPLES) {print $0}}}' \
#  <(awk -F"\t" -v MYPOP=${POP} '{if ($4==MYPOP) print $1}' ${META_DIR}/igsr_samples.tsv) \
#  ${META_DIR}/1kGP.2504_independent_samples.pedigree_info.txt > ${POP}/1KG_${POP}_IND.tsv

cd ${INPUT_DIR}/${POP}

chr=${SLURM_ARRAY_TASK_ID}

PATH_BCFTOOLS="/pasteur/zeus/projets/p02/IGSR/Software/bcftools-1.15"

${PATH_BCFTOOLS}/bcftools view -S 1KG_${POP}_IND.tsv -o CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz ${INPUT_DIR}/ALL/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz

${PATH_BCFTOOLS}/bcftools index CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz -t

# convert using relate software
PATH_RELATE="/pasteur/zeus/projets/p02/IGSR/Software/relate_v1.1.8_x86_64_static/bin"

${PATH_RELATE}/RelateFileFormats --mode ConvertFromVcf --haps CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz.haps --sample CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased_${POP}.vcf.gz.sample -i CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased_${POP}
###############################################################################

