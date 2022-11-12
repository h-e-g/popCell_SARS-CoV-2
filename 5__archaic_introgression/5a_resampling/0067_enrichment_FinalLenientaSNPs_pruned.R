## load libs
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(stringr)
library(tidyr)
library(data.table)
library(tictoc)

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
#Theme set for all plots
theme_set(theme_bw(base_size = 8) +
            theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), strip.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 7), legend.text = element_text(size = 7),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

## wd
# setwd("/Volumes/evo_immuno_pop/users/Javier/")

## vars
temp=commandArgs(TRUE)
num_test=as.numeric(temp[1])
num_resamples=as.numeric(temp[2])
pop=as.character(temp[3])
archaic=as.character(temp[4])
MAF=as.character(temp[5])
DIST=as.logical(temp[6])
LD=as.logical(temp[7])
PRUNE_BGD=as.logical(temp[8])
ADAPTIVE=as.logical(temp[9])
PRUNE_SET=PRUNE_BGD
RESAMP_FROM_SET=TRUE
PRUNE_EUB=FALSE
# num_test <- 100
# num_resamples <- 100
# pop <- "CEU"
# archaic <- "Archaic"
# MAF  <- "Global"
# DIST <- TRUE
# LD <- TRUE

## snps sets
snp_sets <- fread("data/snp_sets/snpSets_June2022.txt") %>% select(set,snps)
# length(unique(snp_sets$set))
# [1] 6
test_name <- unique(snp_sets$set)[num_test]
adj_name=paste0(MAF,'MAF',ifelse(LD,'_LD',''), ifelse(DIST,'_DIST',''),ifelse(PRUNE_BGD,'_prunedBGD',''))

cat("Archaic Lenient",ifelse(ADAPTIVE,'Adpative',''), archaic, ',', pop, 'test nb:',num_test,':',test_name,"- adjusted on", MAF,'MAF,',ifelse(LD,'LDscore,',''), ifelse(DIST,' and DIST',''), ifelse(PRUNE_BGD,' pruning background',''),'\n')

## original mapping file
tot_snps <- fread("results/snps_ids/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_ID_SNPS_REF_ALT.txt") %>% select(ID,snps)

if(PRUNE_BGD){
  if(PRUNE_EUB){
    LD_pruned_EUB=fread(sprintf("%s/single_cell/resources/references/Quach_et_al_2016/LD_pruned_snps_EUB_r2.80.txt",EIP),header=F)$V1
    tot_snps=tot_snps[snps%chin%c(LD_pruned_EUB,snp_sets[set==test_name,snps]),]

  }else{
    LD_pruned_CEU=fread(sprintf("%s/users/Javier/VCF/prune/tagSNPs_%s_maf0.05pct_LD0.8_allCHR.prune.in",EIP,pop),header=F)$V1
    tot_snps=tot_snps[ID%chin%LD_pruned_CEU | snps%chin% snp_sets[set==test_name,snps],]
  }

}
# colnames(tot_snps)[3:4] <- c("REF_ORI","ALT_ORI")
gc()

## read files with factors to account for in this case all SNPs
all_snps <- fread(paste0("results/bins/",pop,"_ALLCHR_popCellSNVs_nomono_maf0.05_hg38_3BINS_4lenientaSNPs_withGlobalMAF.txt.gz")) %>% select(ID,MAF_BIN_0.01,DIST_BIN,LDscore_BIN,MEAN_MAF_BIN_0.01)

## grab only snps from original mapping file
all_snps <- all_snps[all_snps$ID %in% tot_snps$ID,]

## make bin
if (MAF=='Global'){
  all_snps <- all_snps %>% mutate(BIN=MEAN_MAF_BIN_0.01)
}else{
  all_snps <- all_snps %>% mutate(BIN=MAF_BIN_0.01)
}

if (DIST){
  all_snps <- all_snps %>% mutate(BIN=paste0(BIN,"_",DIST_BIN))
}
if (LD){
  all_snps <- all_snps %>% mutate(BIN=paste0(BIN,"_",LDscore_BIN))
}
# all_snps <- all_snps %>% select(ID,BIN)

## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,]
snp_sets <- left_join(snp_sets,tot_snps,by="snps")
eqtls <- unique(snp_sets$ID)
eqtls_snps <- unique(snp_sets$snps)

if(PRUNE_SET){
  popEIP=case_when(pop=='CEU'~'EUB',
                  pop=='CHS'~'ASH',
                  pop=='YRI'~'AFB',
                  TRUE~'.*')
    if(!grepl('quach',tolower(test_name))){
      eqtl_genotypes=fread(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/All_eQTL_and_reQTL_genotypes.tsv.gz",EIP))
    }else{
      eqtl_genotypes=fread(sprintf('%s/single_cell/resources/references/Quach_et_al_2016/All_eQTL_and_reQTL_genotypes_quach.tsv.gz',EIP))
    }
    eqtls=eqtls[eqtls_snps%chin%unique(eqtl_genotypes$ID)]
    eqtls_snps=eqtls_snps[eqtls_snps%chin%unique(eqtl_genotypes$ID)]
    eqtl_genotypes=dcast(eqtl_genotypes,IID~ID,value.var="Number_of_ALT_alelle")
    R2_geno=cor(as.matrix(eqtl_genotypes[grepl(popEIP,IID),mget(eqtls_snps)]))^2
    w.na=which(apply(is.na(R2_geno),1,mean)>.95)
    eqtls=eqtls[-w.na]
    eqtls_snps=eqtls_snps[-w.na]
    eqtl_Clusters=cutree(hclust(as.dist(1-R2_geno[-w.na,-w.na]),method='complete'),h=0.2)
    DT_eQTL=data.table(ID=eqtls,snps=eqtls_snps,eqtl_Clusters)
    set.seed(0)
    DT_eQTL_ind=DT_eQTL[,.SD[sample(1:.N,1)],by=eqtl_Clusters]
    eqtls=DT_eQTL_ind$ID
    eqtls_snps=DT_eQTL_ind$snps
    rm(DT_eQTL_ind,DT_eQTL,eqtl_Clusters,R2_geno,eqtl_genotypes)
    gc()
  }


## read aSNPs file
if (archaic=="Quach"){
  Map_quach=fread(sprintf("%s/single_cell/resources/references/Quach_et_al_2016/Map_imputed_essential_informations.txt.gz",EIP))
  setnames(Map_quach,'snp.name','snps')
  # add Coordinates b38
  cat('merging quach aSNPs with totSNPs')
  Map_quach=merge(Map_quach[aSNP_R2_EUB==TRUE, .(snps,posID_b37=posID,daf_char_EUB,aSNP,aSNP_R2_EUB)],tot_snps,by='snps')[order(ID),]
  #asnps_quach=Map_quach[aSNP_R2_EUB==TRUE,.(ID=ID,HAP_ORIGIN='ALTAI_QUACH')]
  #asnps <- fread(paste0("results/asnps/Final_lenient_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt"))
  # merge(asnps,asnps_quach,by='ID',all.x=T)[,.N,by=.(HAP_ORIGIN.x,MATCH,Population,is.na(HAP_ORIGIN.y))][Population=='CEU',][order(HAP_ORIGIN.x,MATCH,Population),N/sum(N),by=.(HAP_ORIGIN.x,MATCH,Population)]
  # 20 of archaic & Neanderthal SNPs are missing in quach data
  # 40 of AMH snps on archaic & neanderthal haplotypes are missing in quach data
  asnps <- Map_quach[,.(ID=ID,HAP_ORIGIN='ALTAI_QUACH')]
  rm(Map_quach)
}else{
  asnps <- fread(paste0("results/asnps/Final_lenient_Adaptive_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt"))
  if (archaic=="Vindija33.19"){
    archaic2 <- "NEAND"
    }
  if (archaic=="Denisova") {
    archaic2 <- "DENI"
  }
  if (archaic=="Archaic" | archaic=='Any') {
    archaic2 <- toupper(archaic)
  }
  asnps <- asnps %>% filter(Population==pop & HAP_ORIGIN==archaic2)
  if(ADAPTIVE){
    asnps <- asnps %>% filter(Q99==TRUE)
  }
  asnps <- asnps %>% select(ID,HAP_ORIGIN)
  gc()
}
rm(tot_snps,snp_sets)
gc()

## add aSNPS to all_snps
all_snps <- left_join(all_snps,asnps[,c("ID","HAP_ORIGIN")],by="ID")
all_snps$SCORE <- 0
all_snps$HAP_ORIGIN <- ifelse(is.na(all_snps$HAP_ORIGIN),F,T)
#fwrite(all_snps[,table(HAP_ORIGIN,LDscore_BIN,exclude='')],file='results/DistribLDscoreBIN_byArchaicAncestry.tsv',sep='\t')

all_snps <- all_snps %>% select(ID,BIN,SCORE,HAP_ORIGIN)
rm(asnps);gc()

## check if eQTLs are present in file
eqtls_final <- eqtls[(eqtls %in% all_snps$ID)]
rm(eqtls)

## change columns to SCORE and SELECTED
colnames(all_snps) <- c("ID","BIN","SCORE","SELECTED")

## load function
source("scripts/resample_SNPs_matchedbyBIN_SCORE_SELECTED_withDepletion_function.R")

## resample
res <- resample_ids(target_snps = eqtls_final, all_snps_bins = all_snps, num_resamples = num_resamples, rm_target_from_all_snps = ifelse(PRUNE_BGD,RESAMP_FROM_SET,TRUE), use_data.table=TRUE)
res <- as_tibble(res)

## save
write.table(res,
            file=paste0("results/resampling_asnps_pruned/",pop,"_",archaic,"_Final",ifelse(ADAPTIVE,'Adaptive',''),"LenientaSNPs_",test_name,"_adj",adj_name,ifelse(RESAMP_FROM_SET,'','_noResampFromSet'),"_NumResamps",format(num_resamples,scientific=F),".txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")
