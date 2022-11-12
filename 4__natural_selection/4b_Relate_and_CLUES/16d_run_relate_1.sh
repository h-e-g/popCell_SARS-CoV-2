#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 50G
#SBATCH -o "/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/log/%J_relate1.log"
#SBATCH -J rlt
#SBATCH --array=1-22

# sbatch ./16d_run_relate_1.sh CHS
# sbatch ./16d_run_relate_1.sh CEU
# sbatch ./16d_run_relate_1.sh YRI

###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38"
OUT_DIR="${EIP}/project/pop_eQTL/data/4_natural_selection/clues"
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
PATH_MAP="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Genetic_Maps/SHAPEIT_Format/b38"

POP=$1

cd ${OUT_DIR}/relate/${POP}

chr=${SLURM_ARRAY_TASK_ID}

${PATH_RELATE}/bin/Relate \
    --mode All \
    -m 1.25e-8 \
    -N 30000 \
    --haps ${INPUT_DIR}/${POP}/relate_input/${VCF_PREFIX}${chr}.filtered.shapeit2-duohmm-phased_${POP}_input.haps.gz \
    --sample ${INPUT_DIR}/${POP}/relate_input/${VCF_PREFIX}${chr}.filtered.shapeit2-duohmm-phased_${POP}_input.sample.gz \
    --map ${PATH_MAP}/chr${chr}.b38.gmap.gz \
    --seed 1 \
    -o relate_output_1_chr${chr}
