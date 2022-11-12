#!/bin/bash
#SBATCH --qos=normal
#SBATCH -p common
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 8G
#SBATCH -o "/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/log/%A_%a_sample_branch_lengths.log"
#SBATCH -J sbz
#SBATCH --array=1-1782

expresp="eQTL" # if reQTL, length array=1505, else if eQTL length array=12753


RERUN_FILE="20280707_array_timeout_rerun_csh.txt" #1782
#RERUN_FILE="20771531_array_timeout_rerun_yri.txt" #4637


###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38/CHS"
OUT_DIR="${EIP}/project/pop_eQTL/data/4_natural_selection/clues"
POPCELL_DATA="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"

LOG_DIR="${OUT_DIR}/relate/log"

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
PATH_MASK="/pasteur/zeus/projets/p02/IGSR/Genome_Masks/1KG/b38/StrictMask"
PATH_MAP="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Genetic_Maps/SHAPEIT_Format/b38"

NUM=${SLURM_ARRAY_TASK_ID}

#rsid=$(sed "${NUM}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f1)
#chr=$(sed "${NUM}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f2)
#chrpos=$(sed "${NUM}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f3)

NUM_RERUN=$(sed "${NUM}q;d" ${LOG_DIR}/${RERUN_FILE} | cut  -f1)
rsid=$(sed "${NUM_RERUN}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f1)
chr=$(sed "${NUM_RERUN}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f2)
chrpos=$(sed "${NUM_RERUN}q;d" ${OUT_DIR}/${expresp}_metadata_unique.tsv | cut  -f3)

pop=CHS
cutoff=2000

sh ./16i_resample_branch_lengths_per_snp.sh ${rsid} ${chr} ${chrpos} ${pop} ${cutoff} ${expresp}
