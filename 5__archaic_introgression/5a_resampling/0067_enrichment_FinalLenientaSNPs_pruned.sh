#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#SBATCH -o /pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/logs/slurm-%A_%a_LenientaSNPs.out
#SBATCH -J resamp
#SBATCH --mem=8G

## libs
module purge
module load R/4.1.0
module load plink/1.90b6.16
module load tabix
module load samtools/1.10
module load vcftools/0.1.16
module load java/13.0.2
module load shapeit4/4.2.1


################ prune SNPs in LD
# module load plink/1.90b6.16
# module load plink/2.00a3
# WORKDIR='/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/'
#
# for LD in 0.8;
# do
#   for POP in CEU CHS YRI ;
#   do
#     rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_allCHR.prune.in
#     rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD0.2_allCHR.prune.in
#     for CHR in `seq 1 22`;
#      do
#      VCFFILE="${WORKDIR}/VCF/CHS_CEU_YRI_chr${CHR}_popCellSNVs_nomono_hg38.vcf.gz"
#      POPFILE=/pasteur/zeus/projets/p02/IGSR/1KG/1KG_VCF_Files/b38/${POP}/1KG_${POP}_IND.tsv
#        POPFILE_TMP=${SCRATCH}/mrotival/1KG_${POP}_IND_tmp.tsv
#       awk -v pop="$POP" '{$2=$1;print}' $POPFILE > $POPFILE_TMP
#         rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD0.2_${CHR}.*
#         rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.in
#         rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.out
#         rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.nosex
#         rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.log
#         plink --vcf ${VCFFILE} --allow-no-sex --keep ${POPFILE_TMP} --maf 0.05 --indep-pairwise 1000 kb 1 ${LD} --out ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}
#         rm $POPFILE_TMP
#         cat ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.in >> ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_allCHR.prune.in
#       done;
#    done;
#  done;
#awk -v pop="$POP" '{$2=$1;$1=pop; print}' $POPFILE > $POPFILE_TMP

## args
# one of MONO NS aquino, MONO reQTL IAV aquino, MONO eQTL IAV aquino, MONO NS quach, MONO reQTL IAV quach, MONO eQTL IAV quach,
num_test=$SLURM_ARRAY_TASK_ID  # 100, 101, 180, 257-259
num_resamples=$1 # 1000 (faster) or 10000
pop=$2
archaic=$3 #Any OR Quach
MAF=$4 #global (Javier) or local (quach)
DIST=$5 #TRUE (Javier) or FALSE (quach)
LD=$6 #TRUE (Javier) or FALSE (quach)
PRUNE_BGD=$7 # TRUE (quach) or FALSE (Javier)
ADAPTIVE=$8 # TRUE (quach) or FALSE (Javier)

## paths
SS="/pasteur/zeus/projets/p02/evoceania/Javier/Softwares_scripts/"
IGSR="/pasteur/zeus/projets/p02/IGSR/"
REFhg38="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/references/RNA/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
SP="/pasteur/appa/scratch/Javier_tmp/"

## wd
cd "/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier"

echo "Starting time: `date`"


 SNPSET_FILE="data/snp_sets/snpSets_June2022.txt"
 testname=$(awk -v NUM="$num_test" '{if ($3==NUM){print $1; exit}}' $SNPSET_FILE)


if [[ $LD == 'TRUE' ]]; then LD_CHAR='_LD'; else LD_CHAR=''; fi
if [[ $MAF == 'Global' ]]; then MAF_CHAR='GlobalMAF'; else MAF_CHAR='LocalMAF'; fi
if [[ $DIST == 'TRUE' ]]; then DIST_CHAR='_DIST'; else DIST_CHAR=''; fi
if [[ $PRUNE_BGD == 'TRUE' ]]; then PRUNE_BGD_CHAR='_prunedBGD'; else PRUNE_BGD_CHAR='';  fi
if [[ $ADAPTIVE == 'TRUE' ]]; then ADAPTIVE_CHAR='Adaptive'; else ADAPTIVE_CHAR='';  fi

adj_name=${MAF_CHAR}${LD_CHAR}${DIST_CHAR}${PRUNE_BGD_CHAR}
echo ${pop} ${archaic} ${ADAPTIVE_CHAR} ${adj_name} ${num_resamples} 

OUTFILE="results/resampling_asnps_pruned/${pop}_${archaic}_Final${ADAPTIVE_CHAR}LenientaSNPs_${testname}_adj${adj_name}_noResampFromSet_NumResamps${num_resamples}.txt"

## check if outpout exists and, if not, run R
if [ ! -f "${OUTFILE}" ]; then
  Rscript ./scripts/0067_enrichment_FinalLenientaSNPs_pruned.R ${num_test} ${num_resamples} ${pop} ${archaic} ${MAF} ${DIST} ${LD} ${PRUNE_BGD} ${ADAPTIVE}
fi

echo "End time: `date`"
