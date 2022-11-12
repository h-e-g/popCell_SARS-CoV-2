#!/bin/bash
#sbatch --array=1-22 --mem=200000 --qos=geh --partition=geh -J Impute_popCell -o "log/logfile_ImputeEIP_%a.log" ./Imputation_Genotype_AFB_EUB_ASH_ind.sh

source /local/gensoft2/adm/etc/profile.d/modules.sh # source the module loading scripts
module purge # remove modules that may have been loaded by mistake
module load samtools/1.10
module load tabix
module load zlib/1.2.11
#module load vcftools/0.1.16
module load plink/1.90b6.16
module load shapeit4/4.2.1
module load graalvm/ce-java8-20.0.0

module load java/1.8.0
module load beagle/5.1

CHR=${SLURM_ARRAY_TASK_ID}

TYPE="chr${CHR}"
mkdir /pasteur/appa/scratch/public/mrotival/
mkdir /pasteur/appa/scratch/public/mrotival/Imputation
SCRATCH="/pasteur/appa/scratch/public/mrotival/Imputation/${TYPE}"
mkdir $SCRATCH
cd $SCRATCH

EIP_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop"
GENO_DIR="${EIP_DIR}/popCell_data/01_GenotypingData"
WORK_DIR="${EIP_DIR}/single_cell/project/pop_eQTL/Imputation"
OUT_DIR="${EIP_DIR}/popCell_data/02_ImputedGenotypes"

IGSR="/pasteur/zeus/projets/p02/IGSR/"
# shapeit4 comptible maps (not used)
MAP_DIR="${EIP_DIR}/single_cell/resources/references/recombination_maps"
MAP_SHAPEIT_CHR="${MAP_DIR}/genetic_maps.b38/chr${CHR}.b38.gmap"
MAP_PLINK_CHR="${MAP_DIR}/plink.GRCh38.map/plink.chr${CHR}.GRCh38.map"
MAP_PLINK_CHR_UCSC="${MAP_DIR}/plink.GRCh38.map/plink.chr${CHR}.GRCh38.ucsc.map"

# convert the plink format Genotypes to VCF format
# <<<UPDATE THIS>>>
GENOFILE="PopCell_1KG_2977x3723840" # <<<UPDATE THIS>>>
VCF_hg38_FILE="Geno_b38_473Ind_3723840snps"
VCF_hg38="${SCRATCH}/${VCF_hg38_FILE}" # <<<UPDATE THIS>>>

echo PopCell > PopCell_select.fam
echo EvoImmunoPop >> PopCell_select.fam
# /pasteur/appa/scratch/public/mrotival/Imputation/chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nodup.vcf.gz

##### reactivate when rerunning from scratch
plink --bfile ${GENO_DIR}/${GENOFILE} --recode vcf --out ${VCF_hg38}_chr${CHR} --chr ${CHR} --keep-fam PopCell_select.fam

awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' ${VCF_hg38}_chr${CHR}.vcf | grep -ve END | grep -ve SVLEN > ${VCF_hg38}_chr${CHR}.clean.vcf

bgzip ${VCF_hg38}_chr${CHR}.clean.vcf
bcftools index ${VCF_hg38}_chr${CHR}.clean.vcf.gz
tabix -p vcf ${VCF_hg38}_chr${CHR}.clean.vcf.gz

# remove duplicates SNPs from genotyped files
bcftools norm --rm-dup all ${VCF_hg38}_chr${CHR}.clean.vcf.gz -O z -o ${VCF_hg38}_chr${CHR}.nodup.vcf.gz
bcftools index ${VCF_hg38}_chr${CHR}.nodup.vcf.gz
tabix -p vcf ${VCF_hg38}_chr${CHR}.nodup.vcf.gz

##### END reactivate when rerunning from scratch
# copy 1kG ph3 data to scratch, index, and copy the bgzip & tabix index to the /1kG/ folder
# REF_DIR="${EIP_DIR}/Martin/Ressources/1000GFrequencies/RawData"
# for i in `seq 1 22`; do
#   echo $i
# REF_CHR="${REF_DIR}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes"
# LOCAL_REF="${SCRATCH}/genotypes_1kG_phase3.chr${i}"
# cp ${REF_CHR}.vcf.gz ${LOCAL_REF}.vcf.gz
# gunzip ${LOCAL_REF}.vcf.gz
# bgzip ${LOCAL_REF}.vcf
# bcftools index ${LOCAL_REF}.vcf.gz
# tabix -p vcf ${LOCAL_REF}.vcf.gz
# cp ${LOCAL_REF}.vcf.gz ${WORK_DIR}/1kG/
# cp ${LOCAL_REF}.vcf.gz.tbi ${WORK_DIR}/1kG/
# done;

REF_1KG="${IGSR}/1KG/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes2"
LOCAL_REF="${SCRATCH}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"

# gunzip ${LOCAL_REF}.vcf.gz
# bgzip ${LOCAL_REF}.vcf
# bcftools index ${LOCAL_REF}.vcf.gz
# tabix -p vcf ${LOCAL_REF}.vcf.gz

# which version ?
# size of the ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes2.vcf.gz file (214.5Mo according to macos) , matches that of files downloaded from
# http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
# but md5sum is different (size differ by a few bytes 214453750 vs 214453881)

# liftOver 1kG to b38 ?


bcftools norm --rm-dup all ${REF_1KG}.vcf.gz -o ${LOCAL_REF}.nodup.vcf

######################################################
###### before imputation, liftover 1KG to hg38 #######
######################################################

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

bcftools norm --rm-dup all ${LOCAL_REF_GRCh38}.vcf.gz -O z -o ${LOCAL_REF_GRCh38}.nodup.vcf.gz
bcftools index ${LOCAL_REF_GRCh38}.nodup.vcf.gz
tabix -p vcf ${LOCAL_REF_GRCh38}.nodup.vcf.gz

# use Shapeit 4 to phase
mkdir ${WORK_DIR}/log
PHASED_CHR="${VCF_hg38}_chr${CHR}_shapeit4"
shapeit4 --seed 123456 --pbwt-depth 8 --thread 8 --input ${VCF_hg38}_chr${CHR}.nodup.vcf.gz --map ${MAP_SHAPEIT_CHR}.gz --reference ${LOCAL_REF_GRCh38}.vcf.gz --log ${WORK_DIR}/log/phased_${CHR}.log --output ${PHASED_CHR}.vcf.gz --region chr${CHR}
#bgzip ${PHASED_CHR}.vcf

bcftools norm --rm-dup all ${PHASED_CHR}.vcf.gz -O z -o ${PHASED_CHR}.nodup.vcf.gz
bcftools index ${PHASED_CHR}.nodup.vcf.gz
tabix -p vcf ${PHASED_CHR}.nodup.vcf.gz

### move the result to WORK_DIR

mkdir ${OUT_DIR}/phased/
cp ${PHASED_CHR}.nodup.* ${OUT_DIR}/phased/
chmod 775 ${OUT_DIR}/phased/${VCF_hg38_FILE}_chr${CHR}_shapeit4*

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

PHASED_CHR="${VCF_hg38}_chr${CHR}_shapeit4"
IMPUTED_CHR="${VCF_hg38}_chr${CHR}_shapeit4_beagle5"


awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else print $0
      }' ${MAP_PLINK_CHR} > ${MAP_PLINK_CHR_UCSC}

beagle ne=20000 seed=123456 nthreads=8 ref=${LOCAL_REF_GRCh38}.nodup.vcf.gz gt=${PHASED_CHR}.nodup.vcf.gz out=${IMPUTED_CHR} chrom=chr${CHR} map=${MAP_PLINK_CHR_UCSC}

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

bcftools filter -i "DR2>0.9 & AF>0.01 & AF<0.99" ${IMPUTED_CHR}.nodup.vcf.gz -o ${IMPUTED_CHR}_filtered.vcf -O v
bgzip ${IMPUTED_CHR}_filtered.vcf

mkdir ${OUT_DIR}/Imputed/
mkdir ${OUT_DIR}/Imputed/b38/
mv ${IMPUTED_CHR}_filtered.vcf.gz ${OUT_DIR}/Imputed/b38/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz
tabix -p vcf ${OUT_DIR}/Imputed/b38/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz
bcftools index ${OUT_DIR}/Imputed/b38/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz

chmod 775 ${OUT_DIR}/Imputed/b38/${GENOFILE}_chr${CHR}_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct*
rm -r ${SCRATCH}
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
