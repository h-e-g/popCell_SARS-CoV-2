#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 20G
#SBATCH -o "/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/log/%J_popsize.log"
#SBATCH -J psz

# sbatch ./16e_estimate_population_sizes.sh CHS
# sbatch ./16e_estimate_population_sizes.sh CEU
# sbatch ./16e_estimate_population_sizes.sh YRI

###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38"
OUT_DIR="${EIP}/project/pop_eQTL/data/4_natural_selection/clues"
POPCELL_DATA="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
META_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Metadata"

module load plink
module load R
module load samtools
module load tabix
module load vcftools/0.1.16

POP=$1

#cd ${INPUT_DIR}/${POP}/relate_input/
#awk -F"\t" '{if (NR==FNR) {POP_SAMPLES[$1]=$1} else { if ($1 in POP_SAMPLES) {print $1" "$4" "$6}}}' ../1KG_${POP}_IND.tsv ${META_DIR}/igsr_samples.tsv > poplabels.tsv
#################################################################################

# prepare input files
PATH_RELATE="/pasteur/zeus/projets/p02/IGSR/Software/relate_v1.1.8_x86_64_static"
VCF_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05_chr"
PATH_ANCESTOR="/pasteur/zeus/projets/p02/IGSR/Human_Ancestor/GRCh38/homo_sapiens_ancestor_GRCh38"
PATH_MASK="/pasteur/zeus/projets/p02/IGSR/Genome_Masks/1KG/b38/PilotMask"
PATH_MAP="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Genetic_Maps/SHAPEIT_Format/b38"

cd ${OUT_DIR}/relate/${POP}

${PATH_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i relate_output_1 \
              -m 1.25e-8 \
              --seed 1 \
              --poplabels ${INPUT_DIR}/${POP}/relate_input/poplabels.txt \
              -o relate_popsize_${POP} \
              --first_chr 1 \
              --last_chr 22
