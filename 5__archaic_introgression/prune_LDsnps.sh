#!/bin/bash
################################################################################
################################################################################
# File name: pruneLD_snps.sh
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: prune SNPs in LD from 1kG
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=5_00
#SBATCH --output=${LOG_DIR}/5_00__%J.log
#SBATCH -J LDpruning
#SBATCH --mem=8G
################################################################################
################################################################################

# Setup
module load plink/1.90b6.16
module load plink/2.00a3
WORKDIR='/DATA/1KG'

for LD in 0.8;
do
  for POP in CEU CHS YRI ;
  do
    rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_allCHR.prune.in
    for CHR in `seq 1 22`;
     do
     VCFFILE="${WORKDIR}/VCF/CHS_CEU_YRI_chr${CHR}_popCellSNVs_nomono_hg38.vcf.gz"
     POPFILE=${WORKDIR}/1KG_${POP}_IND.tsv
     POPFILE_TMP=${SCRATCH}/mrotival/1KG_${POP}_IND_tmp.tsv
      awk -v pop="$POP" '{$2=$1;print}' $POPFILE > $POPFILE_TMP
        rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.in
        rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.out
        rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.nosex
        rm ${WORKDIR}/VCF/prune/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.log
        plink --vcf ${VCFFILE} --allow-no-sex --keep ${POPFILE_TMP} --maf 0.05 --indep-pairwise 1000 kb 1 ${LD} --out ${WORKDIR}/ld_1kG/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}
        rm $POPFILE_TMP
        cat ${WORKDIR}/ld_1kG/tagSNPs_${POP}_maf0.05pct_LD${LD}_${CHR}.prune.in >> ${WORKDIR}/ld_1kG/tagSNPs_${POP}_maf0.05pct_LD${LD}_allCHR.prune.in
      done;
   done;
 done;
awk -v pop="$POP" '{$2=$1;$1=pop; print}' $POPFILE > $POPFILE_TMP
