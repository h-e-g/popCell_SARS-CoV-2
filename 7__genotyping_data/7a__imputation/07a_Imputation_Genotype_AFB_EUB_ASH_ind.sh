#!/bin/bash
################################################################################
################################################################################
# File name: 07a_Imputation_Genotype_AFB_EUB_ASH_ind.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-chromosome phasing and imputation
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=07a_Impute_popCell
#SBATCH --output=${LOG_DIR}/07a_Impute_popCell__%A_%a.log
#SBATCH --mem=200000
#SBATCH --ntasks=8
#SBATCH --array=1-22
################################################################################
# example
# sbatch ./Imputation_Genotype_AFB_EUB_ASH_ind.sh

source /local/gensoft2/adm/etc/profile.d/modules.sh # source the module loading scripts
module purge # remove modules that may have been loaded by mistake
module load samtools/1.10
module load tabix
module load zlib/1.2.11
module load vcftools/0.1.16
module load plink/1.90b6.16
module load shapeit4/4.2.1
module load graalvm/ce-java8-20.0.0
# load Beagle/5.1
module load java/1.8.0
module load beagle/5.1

CHR=${SLURM_ARRAY_TASK_ID}
DATA_DIR='"../../../DATA'
GENO_DIR="${DATA_DIR}/Genotype/"
WORK_DIR="${GENO_DIR}/Imputation"
TASK_DIR="${WORK_DIR}/chr${CHR}"

OUT_DIR="${GENO_DIR}/Imputed/b38"

mkdir ${GENO_DIR}
mkdir ${GENO_DIR}/Imputation
mkdir ${GENO_DIR}/Imputed
mkdir ${GENO_DIR}/Imputed/b38
mkdir ${TASK_DIR}

cd $TASK_DIR

MAP_DIR="${DATA_DIR}/recombination_maps"
MAP_SHAPEIT_CHR="${MAP_DIR}/genetic_maps.b38/chr${CHR}.b38.gmap"
MAP_PLINK_CHR="${MAP_DIR}/plink.GRCh38.map/plink.chr${CHR}.GRCh38.map"
MAP_PLINK_CHR_UCSC="${MAP_DIR}/plink.GRCh38.map/plink.chr${CHR}.GRCh38.ucsc.map" # USCC Format

# convert the plink format Genotypes to VCF format
GENO_ROOT="Geno_b38_473Ind_3723840snps"
GENOFILE="${GENO_DIR}/${GENO_ROOT}"

echo PopCell > PopCell_select.fam
echo EvoImmunoPop >> PopCell_select.fam

##### reactivate when rerunning from scratch
plink --bfile ${GENOFILE} --recode vcf --out ${GENOFILE}_chr${CHR} --chr ${CHR} --keep-fam PopCell_select.fam

awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' ${GENOFILE}_chr${CHR}.vcf | grep -ve END | grep -ve SVLEN > ${GENOFILE}_chr${CHR}.clean.vcf

bgzip ${GENOFILE}_chr${CHR}.clean.vcf
bcftools index ${GENOFILE}_chr${CHR}.clean.vcf.gz
tabix -p vcf ${GENOFILE}_chr${CHR}.clean.vcf.gz

# remove duplicates SNPs from genotyped files
bcftools norm --rm-dup all ${GENOFILE}_chr${CHR}.clean.vcf.gz -O z -o ${GENOFILE}_chr${CHR}.nodup.vcf.gz
bcftools index ${GENOFILE}_chr${CHR}.nodup.vcf.gz
tabix -p vcf ${GENOFILE}_chr${CHR}.nodup.vcf.gz

##### END reactivate when rerunning from scratch

# 1kG should be downloaded from http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ and added to ${DATA_DIR}/1KG

REF_1KG="${DATA_DIR}/1KG/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"
# gunzip ${REF_1KG}.vcf.gz
# bgzip ${REF_1KG}.vcf
# bcftools index ${REF_1KG}.vcf.gz
# tabix -p vcf ${REF_1KG}.vcf.gz


######################################################
###### before imputation, liftover 1KG to hg38 #######
######################################################

 CHAIN="${DATA_DIR}/liftOver/hg19ToHg38.over.modified.chain"
 GENOME_FASTA="${DATA_DIR}/Genome/RNA/human_iav_sars/fasta/genome"
 GATK="${DATA_DIR}/bin/gatk-4.1.2.0/gatk"
 PICARD="${DATA_DIR}/bin/picard_2.20.1/picard.jar"

java -jar ${PICARD} CreateSequenceDictionary REFERENCE=${GENOME_FASTA}.fa OUTPUT=${GENOME_FASTA}.dict

# remove duplicates
bcftools norm --rm-dup all ${REF_1KG}.vcf.gz -o ${REF_1KG}.nodup.vcf
# decompress, change the chromosome to UCSC format for liftover, compress the original file.
bcftools view ${REF_1KG}.nodup.vcf.gz > ${REF_1KG}.nodup.vcf

awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' ${REF_1KG}.nodup.vcf | grep -ve END | grep -ve SVLEN > ${REF_1KG}.nodup.clean.vcf


REF_1KG_GRCh38="${REF_1KG}.nodup.b38"
REF_1KG_GRCh38_REJECTED="${REF_1KG}.nodup.b38rejected"

# liftover
${GATK} LiftoverVcf --java-options "-Xmx200g" -I=${REF_1KG}.nodup.clean.vcf -O=${REF_1KG_REF_GRCh38}.vcf --CHAIN=${CHAIN} --REJECT=${REF_1KG_GRCh38_REJECTED}.vcf -R=${GENOME_FASTA}.fa --TMP_DIR=$SCRATCH --RECOVER_SWAPPED_REF_ALT=TRUE

bgzip ${REF_1KG_GRCh38}.vcf
bcftools index ${REF_1KG_GRCh38}.vcf.gz
tabix -p vcf ${REF_1KG_GRCh38}.vcf.gz

bcftools norm --rm-dup all ${REF_1KG_GRCh38}.vcf.gz -O z -o ${REF_1KG_GRCh38}.nodup.vcf.gz
bcftools index ${REF_1KG_GRCh38}.nodup.vcf.gz
tabix -p vcf ${REF_1KG_GRCh38}.nodup.vcf.gz

# use Shapeit 4 to phase
mkdir ${WORK_DIR}/log
PHASED_CHR="${TASK_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4"
shapeit4 --seed 123456 --pbwt-depth 8 --thread 8 --input ${GENOFILE}_chr${CHR}.nodup.vcf.gz --map ${MAP_SHAPEIT_CHR}.gz --reference ${REF_1KG_GRCh38}.vcf.gz --log ${WORK_DIR}/log/phased_${CHR}.log --output ${PHASED_CHR}.vcf.gz --region chr${CHR}
bcftools norm --rm-dup all ${PHASED_CHR}.vcf.gz -O z -o ${PHASED_CHR}.nodup.vcf.gz
bcftools index ${PHASED_CHR}.nodup.vcf.gz
tabix -p vcf ${PHASED_CHR}.nodup.vcf.gz

### move the result to WORK_DIR
mkdir ${WORK_DIR}/phased/
cp ${PHASED_CHR}.nodup.* ${WORK_DIR}/phased/
chmod 775 ${WORK_DIR}/phased/${GENO_ROOT}_chr${CHR}_shapeit4*

##### END reactivate when rerunning from scratch

# --thread 8
# time for 400 ind + 2500 1kGph3, no multi-threading
# --pbwt-depth=8 : 1Mb - 37s; 2Mb 172s; 3Mb - 367s
# --pbwt-depth=4 : 1Mb - 26s/44s; 2Mb XXs; 3Mb - 205s ; 4Mb - 235 s; 6Mb - 364 s;  ~ 60s per Mb

# check that phased VCFs ends at the same place as the original
# for i in `seq 1 22`; do
#   tail -n 1 Geno_402Ind_imputed_hg19_chr${i}_UCSCformat.vcf | cut -f1-10
#   tail -n 1 Geno_402Ind_3629019snps_hg19_chr${i}.vcf | cut -f1-10
# done;

# rm ${PHASED_CHR}.vcf.gz

IMPUTED_CHR="${TASK_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4_beagle5"

awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else print $0
      }' ${MAP_PLINK_CHR} > ${MAP_PLINK_CHR_UCSC}

beagle ne=20000 seed=123456 nthreads=8 ref=${REF_1KG_GRCh38}.nodup.vcf.gz gt=${PHASED_CHR}.nodup.vcf.gz out=${IMPUTED_CHR} chrom=chr${CHR} map=${MAP_PLINK_CHR_UCSC}
# time for 400 ind + 2500 1kGph3, no multi-threading
#  1Mb - 40s; 2Mb 50s; 4Mb - 80s; 8Mb 160s; 18 Mb 8 min

# After imputation, filter out low MAF and low quality SNPs
tabix -p vcf ${IMPUTED_CHR}.vcf.gz
bcftools norm --rm-dup all ${IMPUTED_CHR}.vcf.gz -O z -o ${IMPUTED_CHR}.nodup.vcf.gz
bcftools index ${IMPUTED_CHR}.nodup.vcf.gz

bcftools filter -i "DR2>0.9 & AF>0.01 & AF<0.99" ${IMPUTED_CHR}.nodup.vcf.gz -o ${IMPUTED_CHR}_filtered.vcf -O v
bgzip ${IMPUTED_CHR}_filtered.vcf

# move imputed files to ${OUT_DIR}
mv ${IMPUTED_CHR}_filtered.vcf.gz ${OUT_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz
tabix -p vcf ${OUT_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz
bcftools index ${OUT_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz

chmod 775 ${OUT_DIR}/${GENO_ROOT}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct*
