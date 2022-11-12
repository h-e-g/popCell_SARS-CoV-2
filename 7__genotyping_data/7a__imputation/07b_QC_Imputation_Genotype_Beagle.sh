#!/bin/bash
#sbatch --array=1-22 --mem=200000 --qos=geh --partition=geh -J reImpute_popCell -o "log/logfile_reImputeEIP_%a.log" ./QC_Imputation_Genotype_Beagle.sh
#sbatch --array=23-2200 --mem=200000 --qos=geh --partition=geh -J reImpute_popCell -o "log/logfile_reImputeEIP_%a.log" ./QC_Imputation_Genotype_Beagle.sh

source /local/gensoft2/adm/etc/profile.d/modules.sh # source the module loading scripts
module purge # remove modules that may have been loaded by mistake
module load samtools/1.10
module load tabix
module load zlib/1.2.11
module load vcftools/0.1.16
module load plink/1.90b6.16
module load shapeit4/4.2.1
module load graalvm/ce-java8-20.0.0

module load java/1.8.0
module load beagle/5.1


TASK_NB=${SLURM_ARRAY_TASK_ID}

# for TASK_NB in `seq 1 2222`;
# do
SNP_SET=`echo \($TASK_NB-1\)/22 | bc`
CHR=`echo $TASK_NB-22*\($SNP_SET\) | bc`
SNP_SET=`echo $SNP_SET+1 | bc`

TYPE="chr${CHR}_SET${SNP_SET}"
mkdir /pasteur/appa/scratch/public/mrotival/
mkdir /pasteur/appa/scratch/public/mrotival/Imputation
SCRATCH="/pasteur/appa/scratch/public/mrotival/Imputation/${TYPE}"
mkdir $SCRATCH
cd $SCRATCH
SCRATCH_1="/pasteur/appa/scratch/public/mrotival/Imputation/chr${CHR}_SET1"


EIP_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop"
GENO_DIR="${EIP_DIR}/popCell_data/01_GenotypingData"
WORK_DIR="${EIP_DIR}/single_cell/project/pop_eQTL/Imputation"
OUT_DIR="${EIP_DIR}/popCell_data/02_ImputedGenotypes"

IGSR="/pasteur/zeus/projets/p02/IGSR/"
# shapeit4 comptible maps (not used)
MAP_DIR="${EIP_DIR}/single_cell/resources/references/recombination_maps"
MAP_PLINK_CHR_UCSC="${MAP_DIR}/plink.GRCh38.map/plink.chr${CHR}.GRCh38.ucsc.map"

VCF_hg38_FILE="Geno_b38_473Ind_3723840snps"
VCF_hg38="${SCRATCH}/${VCF_hg38_FILE}" # <<<UPDATE THIS>>>

REF_1KG="${IGSR}/1KG/KG_VCF_Files/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes2"
LOCAL_REF="${SCRATCH}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"

PHASED_CHR="${OUT_DIR}/phased/${VCF_hg38_FILE}_chr${CHR}_shapeit4"
PHASED_CHR_PROV="${SCRATCH}/phased_${VCF_hg38_FILE}_chr${CHR}_shapeit4"
IMPUTED_CHR="${VCF_hg38}_chr${CHR}_shapeit4_beagle5"

######################################################
###### before imputation, liftover 1KG to hg38 #######
######################################################

##### extract the SNPs from the desired SNP set
# extract sets of SNP to work with (only those genotyped)
bcftools view -H ${PHASED_CHR}.nodup.vcf.gz  | cut -f 1-5 > ${SCRATCH}/SNP_LIST_FULL.txt
awk '{print $1":"$2}' ${SCRATCH}/SNP_LIST_FULL.txt > ${SCRATCH}/SNP_LIST_FULL_chrpos.txt
sort ${SCRATCH}/SNP_LIST_FULL_chrpos.txt -o ${SCRATCH}/SNP_LIST_FULL_chrpos_srt.txt


# extract sets of SNP to impute from
# delete 1 line from FILE every 100 lines starting at line SNP_SET
sed "${SNP_SET}~100d" ${SCRATCH}/SNP_LIST_FULL.txt > ${SCRATCH}/SNP_LIST_FILTERED.txt
awk '{print $1":"$2}' ${SCRATCH}/SNP_LIST_FILTERED.txt > ${SCRATCH}/SNP_LIST_FILTERED_chrpos.txt
sort ${SCRATCH}/SNP_LIST_FILTERED_chrpos.txt -o ${SCRATCH}/SNP_LIST_FILTERED_chrpos_srt.txt

# extract sets of SNP to filter when reimputing
comm -23 ${SCRATCH}/SNP_LIST_FULL_chrpos_srt.txt ${SCRATCH}/SNP_LIST_FILTERED_chrpos_srt.txt > ${SCRATCH}/SNP_LIST_TO_FILTER_chrpos.txt
comm -23 ${SCRATCH}/SNP_LIST_FULL.txt ${SCRATCH}/SNP_LIST_FILTERED.txt > ${SCRATCH}/SNP_LIST_TO_FILTER.txt

if [ ${TASK_NB} -le 22 ]; then do;
  bcftools norm --rm-dup all ${REF_1KG}.vcf.gz -o ${LOCAL_REF}.nodup.vcf

  GATK="${IGSR}/Automated_Pipeline_For_WGS/Resources/Software/bin/gatk-4.1.2.0/gatk"
  CHAIN="${IGSR}/Chain_Files_For_Liftover/hg19ToHg38/hg19ToHg38.over.modified.chain"
  GENOME_FASTA="${EIP_DIR}/single_cell/resources/references/RNA/human_iav_sars/fasta/genome"
  PICARD="${IGSR}/Automated_Pipeline_For_WGS/Resources/Software/bin/picard_2.20.1/picard.jar"

#java -jar ${PICARD} CreateSequenceDictionary REFERENCE=${GENOME_FASTA}.fa OUTPUT=${GENOME_FASTA}.dict

# decompress, change the chromosome to UCSC format for liftover, compress the original file.
# bcftools view ${LOCAL_REF}.nodup.vcf.gz > ${LOCAL_REF}.nodup.vcf

  awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' ${LOCAL_REF}.nodup.vcf | grep -ve END | grep -ve SVLEN > ${LOCAL_REF}.nodup.clean.vcf

  LOCAL_REF_GRCh38="${LOCAL_REF}.nodup.b38"
  LOCAL_REF_GRCh38_REJECTED="${LOCAL_REF}.nodup.b38rejected"

# liftover
  ${GATK} LiftoverVcf --java-options "-Xmx200g" -I=${LOCAL_REF}.nodup.clean.vcf -O=${LOCAL_REF_GRCh38}.vcf --CHAIN=${CHAIN} --REJECT=${LOCAL_REF_GRCh38_REJECTED}.vcf -R=${GENOME_FASTA}.fa --TMP_DIR=$SCRATCH --RECOVER_SWAPPED_REF_ALT=TRUE
# current ERROR:
# java.lang.IllegalStateException: Key SVLEN found in VariantContext field INFO at chr22:16437285 but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.

  bgzip ${LOCAL_REF_GRCh38}.vcf
  bcftools index ${LOCAL_REF_GRCh38}.vcf.gz
  tabix -p vcf ${LOCAL_REF_GRCh38}.vcf.gz

  bcftools view -T ${SCRATCH}/SNP_LIST_FULL.txt ${LOCAL_REF_GRCh38}.vcf.gz -O z -o ${LOCAL_REF_GRCh38}_filtered.vcf.gz
  bcftools norm --rm-dup all ${LOCAL_REF_GRCh38}_filtered.vcf.gz -O z -o ${LOCAL_REF_GRCh38}.nodup.vcf.gz
  bcftools index ${LOCAL_REF_GRCh38}.nodup.vcf.gz
  tabix -p vcf ${LOCAL_REF_GRCh38}.nodup.vcf.gz

  done;
else do;
  LOCAL_REF_1="${SCRATCH_1}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"
  cp ${LOCAL_REF_1}.nodup.vcf.gz ${SCRATCH}
fi;

bcftools view -T ${SCRATCH}/SNP_LIST_FILTERED.txt ${PHASED_CHR}.nodup.vcf.gz -O z -o ${PHASED_CHR_PROV}.nodup_filtered.vcf.gz
bcftools index ${PHASED_CHR_PROV}.nodup_filtered.vcf.gz
tabix -p vcf ${PHASED_CHR_PROV}.nodup_filtered.vcf.gz

beagle ne=20000 seed=123456 nthreads=8 ref=${LOCAL_REF_GRCh38}.nodup.vcf.gz gt=${PHASED_CHR_PROV}.nodup_filtered.vcf.gz out=${IMPUTED_CHR} chrom=chr${CHR}=${START}-${END} map=${MAP_PLINK_CHR_UCSC}
# excludemarkers=${SCRATCH}/SNP_LIST_TO_FILTER_chrpos.txt
# ExcludeList="${WORK_DIR}/Held_out_SNPs_chr1.txt"

# run Beagle5

# time for 400 ind + 2500 1kGph3, no multi-threading
#  1Mb - 40s; 2Mb 50s; 4Mb - 80s; 8Mb 160s; 18 Mb 8 min

##############################################################################
###                QC: Evaluate the quality of the imputation              ###
###       based on a set of 500 SNPs that held out for imputation          ###
##############################################################################

#### perform imputation holding SNPs out
# beagle ne=20000 seed=123456 nthreads=1 ref=${LOCAL_REF}.vcf.gz gt=${WORK_DIR}/All_but_Held_out_SNPs_chr1.vcf out=${IMPUTED_CHR}_Held_out_SNPs chrom=${CHR} map=${MAP_PLINK_CHR} ap=true gp=true

# intersect all imputed SNPs with the set of genotyped SNPs that we held out.
#bgzip ${WORK_DIR}/Held_out_SNPs_chr1.vcf
#bcftools sort ${WORK_DIR}/Held_out_SNPs_chr1.vcf.gz -o ${WORK_DIR}/Held_out_SNPs_chr1.vcf.gz -Oz
#bcftools index ${WORK_DIR}/Held_out_SNPs_chr1.vcf.gz
#bcftools index ${IMPUTED_CHR}_Held_out_SNPs.vcf.gz

#bcftools isec -p ${SCRATCH}/ImputationQC -Oz ${WORK_DIR}/Held_out_SNPs_chr1.vcf.gz ${IMPUTED_CHR}_Held_out_SNPs.vcf.gz

#### the resulting intersect will be used for assessing the quality of the imputation
#### (ground truth VS imputation)

##############################################################################
###                                 END QC                                 ###
##############################################################################
# tabix -p vcf ${IMPUTED_CHR}.vcf.gz
# rm ${LOCAL_REF}.vcf.gz

# After imputation, filter out low MAF and low quality SNPs

bcftools norm --rm-dup all ${IMPUTED_CHR}.vcf.gz -O z -o ${IMPUTED_CHR}.nodup.vcf.gz
bcftools index ${IMPUTED_CHR}.nodup.vcf.gz



bcftools view -T ${SCRATCH}/SNP_LIST_TO_FILTER.txt ${IMPUTED_CHR}.nodup.vcf.gz -O z -o ${IMPUTED_CHR}.nodup_snpSet${SNP_SET}.vcf.gz

GENOFILE=""
mkdir ${OUT_DIR}/Imputed/b38/reImputationQC/
mv ${IMPUTED_CHR}.nodup_filtered_${SNP_SET}.vcf.gz ${OUT_DIR}/Imputed/b38/reImputationQC/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_snpSet${SNP_SET}.vcf.gz
tabix -p vcf ${OUT_DIR}/Imputed/b38/reImputationQC/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_snpSet${SNP_SET}.vcf.gz
bcftools index ${OUT_DIR}/Imputed/b38/reImputationQC/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_snpSet${SNP_SET}.vcf.gz

#bcftools filter -i "DR2>0.9 & AF>0.01 & AF<0.99" ${IMPUTED_CHR}.nodup.vcf.gz -o ${IMPUTED_CHR}_filtered.vcf -O v
#bgzip ${IMPUTED_CHR}_filtered.vcf

chmod 775 ${OUT_DIR}/Imputed/b38/reImputationQC/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_snpSet${SNP_SET}*
if ${TASK_NB}>22;
do;
  rm -r ${SCRATCH}
done;
########################################################################
#############################   OBSOLETE    ############################
###### After imputation + filtering, liftover everything to hg38 #######
########################################################################

# GATK="/pasteur/entites/Geh/Shared_Programs/WES_WGS/SOFTWARE/bin/gatk-4.1.2.0/gatk"
# CHAIN="/pasteur/entites/Geh/Shared_Programs/WES_WGS/DATA/5_Data_for_Liftover/hg19ToHg38.over.modified.chain"
# GENOME_FASTA="/pasteur/projets/policy01/evo_immuno_pop/single_cell/resources/references/RNA/human_iav_sars/fasta/genome"
# PICARD="/pasteur/entites/Geh/Shared_Programs/WES_WGS/SOFTWARE/bin/picard/2.20.1/picard.jar"

# java -jar ${PICARD} CreateSequenceDictionary REFERENCE=${GENOME_FASTA}.fa OUTPUT=${GENOME_FASTA}.dict

# decompress, change the chromosome to UCSC format for liftover, compress the original file.
# mv ${IMPUTED_CHR}_filtered.vcf.gz ${IMPUTED_CHR}_filtered.vcf
# awk '{
#         if($0 !~ /^#/)
#             print "chr"$0;
#         else if(match($0,/(##contig=<ID=)(.*)/,m))
#             print m[1]"chr"m[2];
#         else print $0
#       }' ${IMPUTED_CHR}_filtered.vcf > ${IMPUTED_CHR}_UCSCformat.vcf
#
# grep -ve END ${IMPUTED_CHR}_UCSCformat.vcf > ${IMPUTED_CHR}_clean.vcf
#
# VCF_GRCh38="${SCRATCH}/Geno_402Ind_imputed_liftedhg38_chr${CHR}"
# VCF_GRCh38rejected="${SCRATCH}/Geno_402Ind_imputed_rejected_chr${CHR}"
#
# # liftover
# ${GATK} LiftoverVcf --java-options "-Xmx200g" -I=${IMPUTED_CHR}_clean.vcf -O=${VCF_GRCh38}.vcf --CHAIN=${CHAIN} --REJECT=${VCF_GRCh38rejected}.vcf -R=${GENOME_FASTA}.fa --TMP_DIR=$SCRATCH --RECOVER_SWAPPED_REF_ALT=TRUE

# compress b37 Imputed file and move to EVO_IMMUNO_POP
# bgzip ${IMPUTED_CHR}_filtered.vcf
# mkdir ${WORK_DIR}/b37/
# mv ${IMPUTED_CHR}_filtered.vcf.gz ${WORK_DIR}/b37/

# compress b38 Imputed file and move to EVO_IMMUNO_POP

###### How to remove chr* prefixes from seqnames
# zcat ${VCF_hg19}_chr.vcf.gz | awk '{gsub(/^chr/,""); print}' | bgzip -c > ${VCF_hg19}_NoChrTMP.vcf.gz
# tabix -p vcf ${VCF_hg19}_NoChrTMP.vcf.gz
#
# ${GATK} --java-options "-Xmx8G -XX:ParallelGCThreads=1" UpdateVCFSequenceDictionary \
#       -V  ${VCF_hg19}_NoChrTMP.vcf.gz \
#       --output  ${VCF_hg19}_NoChr.vcf.gz \
#       --source-dictionary ${GENOME_DICT}.dict \
#       --replace=true

####################################################################################
######## Failed attempt at using Impute5 on tars/ replaced by beagle5.1 ############
####################################################################################
# convert to impute2 format
# bcftools convert --haplegendsample ${PHASED_CHR} --vcf-ids ${PHASED_CHR}.vcf.gz
# bcftools convert --haplegendsample ${LOCAL_REF} --vcf-ids ${LOCAL_REF}.vcf.gz

# docker pull ubuntu
# run the docker
# docker run -it --rm --name impute5 -v  ubuntu
# singularity pull docker://ubuntu
# singularity run --home $HOME:/home ubuntu_latest.sif

#run the docker
#[sudo] docker run -it --rm --name impute5 -v <impute_folder>:/home/impute5 ubuntu

#from the new command line, simply call impute5:
#home/impute5/bin/impute5_v1.1.3_static --h /home/impute5/test/reference.bcf --g /home/impute5/test/target.bcf --m /home/impute5/test/chr20.b37.gmap.gz --o home/impute5/test/impute.bcf --r 20:1000000-4000000

# docker run -it --rm --name impute -v /Users/mrotival/WORK/02_data/Imputation:/home ubuntu


#from the new command line, simply call impute5:
#  for i in `seq 2 22`; do
#  #mv $WORK_DIR/b37/Geno_402Ind_imputed_hg19_chr${i}.vcf.gz $WORK_DIR/b37/Geno_402Ind_phased_hg19_chr${i}.vcf.gz
#  tabix -p vcf $WORK_DIR/b37/Geno_402Ind_phased_hg19_chr${i}.vcf.gz
#  done;
#
# /home/impute5/impute5_v1.1.3/impute5_v1.1.3_static --h /home/1kG/genotypes_1kG_phase3.chr1.vcf.gz --g /home/b37/Geno_402Ind_phased_hg19_chr1.vcf.gz --m /home/genetic_maps.b37_SHAPEIT4/chr1.b37.gmap.gz --o /home/b37/Geno_402Ind_imputed_hg19_chr1.vcf.gz --r 1 --b 1000
#
# for i in `seq 2 22`; do
#   /home/impute5/impute5_v1.1.3/impute5_v1.1.3_static --h /home/1kG/genotypes_1kG_phase3.chr${i}.vcf.gz --g /home/b37/Geno_402Ind_phased_hg19_chr${i}.vcf.gz --m /home/genetic_maps.b37_SHAPEIT4/chr${i}.b37.gmap.gz --o /home/b37/Geno_402Ind_imputed_hg19_chr${i}.vcf.gz --r ${i}  --b 1000
# done;

# run Impute5
# ${WORK_DIR}/impute5_v1.1.3/impute5_v1.1.3_static --seed 123456 --h ${LOCAL_REF}.vcf.gz --g ${PHASED_CHR}.vcf.gz --o ${IMPUTED_CHR}.bcf --r 1:2000000-3000000 --m ${MAP_CHR}.gz
####################################################################################
################ END of failed attempt at using Impute5 on tars      ###############
# ####################################################################################
# for CHR in `seq 1 22`;
#   do
#     VCF_b37="${WORK_DIR}/b37/Geno_402Ind_imputed_hg19_chr${CHR}"
#     IMPUTED_CHR="${SCRATCH}/Geno_402Ind_imputed_hg19_chr${CHR}"
# #    mv ${IMPUTED_CHR}.vcf.gz ${VCF_b37}.vcf.gz
# #    mv ${IMPUTED_CHR}.vcf.gz.csi ${VCF_b37}.vcf.gz.csi
#     mv ${IMPUTED_CHR}_filtered.vcf.gz ${VCF_b37}_filtered.vcf.gz
#   done;
