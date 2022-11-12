#!/bin/bash
#SBATCH --qos=geh
#SBATCH -p gehbigmem
#SBATCH --mail-user=yaaquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --mem 50G
#SBATCH -o "/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/log/%A_%a_sample_branch_lengths.log"
#SBATCH -J sbz

###############################################################################
                             # Prep files #
###############################################################################

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell"
INPUT_DIR="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38/CHS"
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
PATH_CLUES="/pasteur/zeus/projets/p02/IGSR/Software/clues"
VCF_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05_chr"
PATH_ANCESTOR="/pasteur/zeus/projets/p02/IGSR/Human_Ancestor/GRCh38/homo_sapiens_ancestor_GRCh38"
PATH_MASK="/pasteur/zeus/projets/p02/IGSR/Genome_Masks/1KG/b38/StrictMask"
PATH_MAP="/pasteur/zeus/projets/p02/IGSR/1KG/1KG_Genetic_Maps/SHAPEIT_Format/b38"


rsid=$1
chr=$2
chrpos=$3
POP=$4
CUTOFF=$5
EXPRESP=$6

echo "Re-sampling branch lengths. ${chr}:${chrpos} (${rsid}) in ${POP} (${EXPRESP}, CutOff=${CUTOFF})."

if test -f "${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/${EXPRESP}/${rsid}/relate_samplebranch_${rsid}.timeb"; then
    echo "relate_samplebranch_${rsid}.timeb  exists."
else
  cd ${OUT_DIR}/relate/${POP}
  mkdir -p ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/${EXPRESP}/${rsid}

  ${PATH_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
          -i relate_branchlengths_updated5_${POP}_${chr} \
          -o ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/${EXPRESP}/${rsid}/relate_samplebranch_${rsid} \
          -m 1.25e-8 \
          #--coal spiedel_popsize_${POP}.coal \
          --coal relate_popsize_updated4_${POP}.coal
          --format b \
          --first_bp ${chrpos} \
          --last_bp ${chrpos} \
          --num_samples 100 
fi

cd ${OUT_DIR}/${POP}_CutOff${CUTOFF}/SNPs/${EXPRESP}/${rsid}

echo "Running CLUES. ${chr}:${chrpos} (${rsid})."

source /pasteur/zeus/projets/p02/IGSR/Software/clues_conda/bin/activate

python3 ${PATH_CLUES}/inference_modified_path.py --times relate_output_2_${rsid} --coal ${OUT_DIR}/relate/${POP}/relate_popsize_updated4_${POP}.coal --out clues_output_2_${rsid} --tCutoff ${CUTOFF}



################################################################################
# Example of CHS
# about half the jobs failed. trying to find out why.
#LOG_PATH="${EIP}/project/pop_eQTL/data/4_natural_selection/clues/relate/log"

#sacct -j 18002624 -o JobID%-20 | grep 18002624 | grep -v batch | grep -v ext | wc -l
#12753

#sacct -j 18002624 -o JobID%-20,Partition,State | \
#  grep FAILED | grep common | \
#  cut -f1 -d" " > ${LOG_PATH}/18002624_array_failed_taskids.txt
#
#wc -l ${LOG_PATH}/18002624_array_failed_taskids.txt
#4598 ../../relate/log/18002624_array_failed_taskids.txt

#cat /pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/log/18002624*_sample_branch_lengths.log | grep "Error: No SNPs are mapping to tree" | wc -l
#4598

# 4598/12753 jobs failed
# 4589 jobs failed because the SNP was not found in the relate (round 2) output

# which SNPs?

#array=$(cat ${LOG_PATH}/18002624_array_failed_taskids.txt | sed 's/18002624_//g')
#
#for j in ${array[@]};
#  do
#    grep -e "^Running CLUES." ${LOG_PATH}/18002624_${j}_sample_branch_lengths.log | sed "s/Running CLUES\. //g" | sed "s/\.$//g" | sed "s/(//g" | sed "s/)//g" | sed "s/ /\t/g"
#  done > ${LOG_PATH}/18002624_failed_SNPs.txt

# I merged the SNP rsid with the eQTL info
#${LOG_PATH}/18002624_failed_SNPs_info.tsv.gz
# 3132 unique genes
#> failed_clues_info[,table(celltype,state)]
#                    state
#celltype             COV IAV  NS
#  B                  339 364 339
#  B.M.K              239 233 220
#  B.M.L               29  55  48
#  B.N.K              129 138 120
#  B.N.L               46 104  80
#  cDC                 19  15  12
#  ILC                  3   1   6
#  MAIT               131 137 183
#  MONO               627 301 689
#  MONO.CD14          607 255 658
#  MONO.CD14.INFECTED   0 107   0
#  MONO.CD16           80  92 115
#  NK                 268 357 297
#  NK.CD56brt          16  56   7
#  NK.CD56dim         194 200 217
#  NK.M.LIKE           64 117  30
#  pDC                 14   7  12
#  Plasmablast          3   1  10
#  T.CD4              751 850 795
#  T.CD4.E            465 580 558
#  T.CD4.N            527 596 520
#  T.CD8              620 661 657
#  T.CD8.CM.EM        256 309 264
#  T.CD8.EMRA         179 159 148
#  T.CD8.N            370 368 339
#  T.gd                 5  11   7
#  T.Reg               72 149  99

# Are this SNPs in the genomic mask?
# 4057/4597 are in the genomic mask

# What about the 540 that are not in the mask?

#> failed_clues_info[chrpos%in%failed_clues_not_in_mask_vec,length(unique(gene_name))]
#[1] 462
#  B                   44  52  51
#  B.M.K               36  41  35
#  B.M.L                4   8   7
#  B.N.K               17  18  12
#  B.N.L               10  10  15
#  cDC                  7   7   2
#  ILC                  1   0   2
#  MAIT                17  22  36
#  MONO                66  30  80
#  MONO.CD14           62  34  78
#  MONO.CD14.INFECTED   0   9   0
#  MONO.CD16           12   6  19
#  NK                  35  54  37
#  NK.CD56brt           4  12   2
#  NK.CD56dim          21  29  41
#  NK.M.LIKE           13  22   7
#  pDC                  5   3   4
#  Plasmablast          1   0   3
#  T.CD4              103 108 104
#  T.CD4.E             60  76  81
#  T.CD4.N             68  78  82
#  T.CD8               93  96  95
#  T.CD8.CM.EM         36  44  36
#  T.CD8.EMRA          24  27  30
#  T.CD8.N             44  42  45
#  T.gd                 1   1   1
#  T.Reg               14  22  16

