################################################################################
################################################################################
# File name: 5a1_enrichment_introgressedSNPs_pruned.R
# Author:  J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute excess of introgressed SNPs among various eQTL lists
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

## snps sets
snp_sets <- fread(sprintf("%s/All_eQTL_snpsSets.txt.gz",EQTL_DIR)) %>% select(set,snps)
test_name <- unique(snp_sets$set)[num_test]

cat("Introgressed SNP, ", pop, 'test nb:',num_test,':',test_name,"/n")

## original mapping file
SNP_info=getMap(annotate=FALSE,ARCHAIC=TRUE,popgen=TRUE)

#define MAF bins
SNP_info[,MAF_YRI:=pmin(DAF_or_MAF_YRI,1-DAF_or_MAF_YRI)]
SNP_info[,MAF_CEU:=pmin(DAF_or_MAF_CEU,1-DAF_or_MAF_CEU)]
SNP_info[,MAF_CHS:=pmin(DAF_or_MAF_CHS,1-DAF_or_MAF_CHS)]
SNP_info[,MAF_SNP:=case_when(pop=='CEU'~MAF_CEU,
                              pop=='CHS'~MAF_CHS,
                              TRUE~NA)]
SNP_info[,MAF_BIN_0.01:=cut(MAF_SNP,breaks = seq(0,0.5,by=0.01),include.lowest = T)]

#define Distance bins
SNP_info[,DIST_BIN:=cut(DistGene,breaks = c(0,1000,5000,10000,20000,50000,100000,Inf),include.lowest = T)]

#define LD score bins
if(pop=='CEU'){
  all_snps=SNP_info[MAF_CEU>=0.05,]
  all_snps[,LD_score_BIN:=cut(LD_score_CEU,breaks = quantile(LD_score_CEU, seq(0,1,0.1)),include.lowest = T)]
  all_snps[,SCORE:=0]
  all_snps[,SELECTED:=(HAP_ORIGIN_CEU!='')]
}
if(pop=='CHS'){
  all_snps=SNP_info[MAF_CHS>=0.05,]
  all_snps[,LD_score_BIN:=cut(LD_score_CHS,breaks = quantile(LD_score_CHS, seq(0,1,0.1)),include.lowest = T)]
  all_snps[,SCORE:=0]
  all_snps[,SELECTED:=(HAP_ORIGIN_CHS!='')]
}

## make bin variable
all_snps <- all_snps %>% mutate(BIN=paste0(MAF_BIN_0.01,"_",DIST_BIN,"_",LDscore_BIN))

## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,.(set,rsID=snps)]
snp_sets <- merge(snp_sets,all_snps[,.(rsID,ID)],by="rsID")
eqtls <- unique(snp_sets$ID)
eqtls_snps <- unique(snp_sets$rsID)

############# eQTL pruning  #############
# load and format eQTL genotypes
eqtl_genotypes=fread(sprintf("%s/All_eQTL_and_reQTL_genotypes.tsv.gz",EQTL_DIR))
eqtl_genotypes=dcast(eqtl_genotypes,IID~ID,value.var="Number_of_ALT_alelle")
# check: discard if genotype is missing
eqtls=eqtls[eqtls_snps%chin%unique(eqtl_genotypes$ID)]
eqtls_snps=eqtls_snps[eqtls_snps%chin%unique(eqtl_genotypes$ID)]
# compute r2 between eQTL in the target population
popEIP=case_when(pop=='CEU'~'EUB',
                  pop=='CHS'~'ASH',
                  pop=='YRI'~'AFB',
                  TRUE~'.*')
R2_geno=cor(as.matrix(eqtl_genotypes[grepl(popEIP,IID),mget(eqtls_snps)]))^2
# check : discard eQTL snps if r2 is missing
w.na=which(apply(is.na(R2_geno),1,mean)>.95)
eqtls=eqtls[-w.na]
eqtls_snps=eqtls_snps[-w.na]
# group tighly linked snps
eqtl_Clusters=cutree(hclust(as.dist(1-R2_geno[-w.na,-w.na]),method='complete'),h=0.2)
DT_eQTL=data.table(ID=eqtls,snps=eqtls_snps,eqtl_Clusters)
# keep a single eQTL per group
set.seed(0)
DT_eQTL_ind=DT_eQTL[,.SD[sample(1:.N,1)],by=eqtl_Clusters]
# update eQTL SNP lists
eqtls=DT_eQTL_ind$ID
eqtls_snps=DT_eQTL_ind$snps
rm(DT_eQTL_ind,DT_eQTL,eqtl_Clusters,R2_geno,eqtl_genotypes)
gc()
############# end: eQTL pruning #############

## check if eQTLs are present in file
eqtls_final <- eqtls[(eqtls %in% all_snps$ID)]
rm(eqtls)
gc()

## prune background SNPs
if(pop=='CEU'){
  all_snps <- all_snps[ID %chin% LD_pruned_CEU | rsID%chin% eqtls_snps,]
}
if(pop=='CHS'){
  all_snps <- all_snps[ID %chin% LD_pruned_CHS | rsID%chin% eqtls_snps,]
}

# extract required columns
all_snps=all_snps[,.("ID","BIN","SCORE","SELECTED")]

## load function
source("scripts/resample_SNPs_matchedbyBIN_SCORE_SELECTED_withDepletion_function.R")

## resample
res <- resample_ids(target_snps = eqtls_final, all_snps_bins = all_snps, num_resamples = num_resamples, rm_target_from_all_snps = T, use_data.table=TRUE)
res <- as.data.table(res)

## save
OUT_FILE=sprintf("%s/%s_introgressed_%s_NumTest%s_NumResamp%s.txt",DAT_RESAMP_ASNP_DIR,pop,test_name,num_test,format(num_resamples,scientific=F))
fwrite(res,file=OUT_FILE,sep = "\t")
