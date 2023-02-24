################################################################################
################################################################################
# File name: 4a1_compute_enrichment_PBS_3factors.R
# Author:  J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute excess of high PBS snps among various eQTL lists
# relative to random matched SNPs, and provide a p-value
# Effector script
################################################################################
################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

## load resampling function
source(sprintf("%s/resample_SNPs_matchedbyBIN_SCORE_SELECTED_withDepletion_function.R",MISC_DIR))

## vars
temp=commandArgs(TRUE)
num_test=as.numeric(temp[1])
num_resamples=as.numeric(temp[2])
pop=as.character(temp[3])

# example paramaters for tests
# num_test <- 24
# num_resamples <- 20
# pop <- "CHS"

## snps sets
snp_sets <- fread(sprintf("%s/All_eQTL_snpsSets.txt.gz",EQTL_DIR)) %>% select(set,snps)
# length(unique(snp_sets$set))

test_name <- unique(snp_sets$set)[num_test]

cat("PBS", pop, 'test nb:',num_test,':',test_name,'\n')

## original mapping file
SNP_info=getMap(annotate=FALSE,popgen=TRUE)

#define MAF bins
SNP_info[,MAF_YRI:=pmin(DAF_or_MAF_YRI,1-DAF_or_MAF_YRI)]
SNP_info[,MAF_CEU:=pmin(DAF_or_MAF_CEU,1-DAF_or_MAF_CEU)]
SNP_info[,MAF_CHS:=pmin(DAF_or_MAF_CHS,1-DAF_or_MAF_CHS)]
SNP_info[,MEAN_MAF:=(MAF_YRI+MAF_CEU+MAF_CHS)/3]
SNP_info[,MEAN_MAF_BIN_0.01:=cut(MEAN_MAF,breaks = seq(0,0.5,by=0.01),include.lowest = T)]

#define Distance bins
SNP_info[,DIST_BIN:=cut(DistGene,breaks = c(0,1000,5000,10000,20000,50000,100000,Inf),include.lowest = T)]

#define LD score bins
if(pop=='YRI'){
  all_snps=SNP_info[MAF_YRI>=0.05,]
  all_snps[,LD_score_BIN:=cut(LD_score_YRI,breaks = quantile(LD_score_YRI, seq(0,1,0.1)),include.lowest = T)]
  all_snps[,SCORE:=PBS_YRI]
  all_snps[,SELECTED:=P_PBS_YRI<.01]
}
if(pop=='CEU'){
  all_snps=SNP_info[MAF_CEU>=0.05,]
  all_snps[,LD_score_BIN:=cut(LD_score_CEU,breaks = quantile(LD_score_CEU, seq(0,1,0.1)),include.lowest = T)]
  all_snps[,SCORE:=PBS_CEU]
  all_snps[,SELECTED:=P_PBS_CEU<.01]
}
if(pop=='CHS'){
  all_snps=SNP_info[MAF_CHS>=0.05,]
  all_snps[,LD_score_BIN:=cut(LD_score_CHS,breaks = quantile(LD_score_CHS, seq(0,1,0.1)),include.lowest = T)]
  all_snps[,SCORE:=PBS_CHS]
  all_snps[,SELECTED:=P_PBS_CHS<.01]
}

## make bin
all_snps[,BIN:=paste0(LD_score_BIN,"_",MEAN_MAF_BIN_0.01,"_",DIST_BIN)]

## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,.(set,rsID=snps)]
snp_sets <- merge(snp_sets,all_snps[,.(rsID,ID)],by="rsID")
eqtls <- unique(snp_sets$ID)
rm(snp_sets)
gc()

## check if eQTLs are present in file
eqtls_final <- eqtls[(eqtls %in% all_snps$ID)]
rm(eqtls)
gc()

all_snps=all_snps[,.("ID","BIN","SCORE","SELECTED")]

## resample
res <- resample_ids(target_snps = eqtls_final, all_snps_bins = all_snps, num_resamples = num_resamples, rm_target_from_all_snps = T, use_data.table=TRUE)
res <- as.data.table(res)

## save
OUT_FILE=sprintf("%s/%s_PBS_%s_NumTest%s_NumResamp%s.txt",DAT_RESAMP_DIR,pop,test_name,num_test,format(num_resamples,scientific=F))
fwrite(res,file=OUT_FILE,sep = "\t")
