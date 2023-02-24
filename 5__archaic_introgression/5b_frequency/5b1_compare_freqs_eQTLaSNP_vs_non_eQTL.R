#!/bin/bash
################################################################################
################################################################################
# File name: 5b1_compare_freqs_eQTLaSNP_vs_non_eQTL.R
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compare frequency of (pruned) aSNPs between eQTL and non eQTL
# effector script
################################################################################

# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare usefule functions
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))

## vars
temp=commandArgs(TRUE)
num_test=as.numeric(temp[1])
pop=as.character(temp[2])

## snps sets
snp_sets <- fread(sprintf("%s/All_eQTL_snpsSets.txt.gz",EQTL_DIR)) %>% select(set,snps)
test_name <- unique(snp_sets$set)[num_test]

## read files with factors to account for in this case all SNPs

## original mapping file
SNP_info=getMap(annotate=FALSE,ARCHAIC=TRUE,popgen=TRUE)

#define MAF bins
SNP_info[,MAF_YRI:=pmin(DAF_or_MAF_YRI,1-DAF_or_MAF_YRI)]
SNP_info[,MAF_CEU:=pmin(DAF_or_MAF_CEU,1-DAF_or_MAF_CEU)]
SNP_info[,MAF_CHS:=pmin(DAF_or_MAF_CHS,1-DAF_or_MAF_CHS)]
SNP_info[,MAF_SNP:=case_when(pop=='CEU'~MAF_CEU,
                              pop=='CHS'~MAF_CHS,
                              TRUE~NA)]
SNP_info[,MEAN_MAF:=(MAF_YRI+MAF_CHS+MAF_CEU)/3]

SNP_info[,MAF_BIN_0.01:=cut(MAF_SNP,breaks = seq(0,0.5,by=0.01),include.lowest = T)]


#define LD score bins
if(pop=='CEU'){
  all_snps=SNP_info[MAF_CEU>=0.05,]
  all_snps[,Introgressed:=(HAP_ORIGIN_CEU!='')]
}
if(pop=='CHS'){
  all_snps=SNP_info[MAF_CHS>=0.05,]
  all_snps[,Introgressed:=(HAP_ORIGIN_CHS!='')]
}

## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,.(set,rsID=snps)]
snp_sets <- merge(snp_sets,all_snps[,.(rsID,posID)],by="rsID")
eqtls <- unique(snp_sets$posID)
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
eqtls_final <- eqtls[(eqtls %in% all_snps$posID)]
rm(eqtls)
gc()


## prune background SNPs
if(pop=='CEU'){
  all_snps <- all_snps[,tag_SNP:=posID %chin% LD_pruned_CEU]
}
if(pop=='CHS'){
  all_snps <- all_snps[,tag_SNP:=posID %chin% LD_pruned_CHS]
  }
all_snps[,eQTL:=posID%in%eqtls_final]

archaic_IDs=all_snps[Introgressed==TRUE,unique(posID)]

########################################################################################
########   obtain frequency of introgressed allele the pop under consideration ########
########################################################################################

if(pop=='CEU'){
  aSNP_good=SNP_info[!is.na(ASNP_FREQ_CEU),/(rsID,posID,ORIGIN,ASNP_FREQ=ASNP_FREQ_CEU)]
  }
if(pop=='CHS'){
  aSNP_good=SNP_info[!is.na(ASNP_FREQ_CHS),/(rsID,posID,ORIGIN,ASNP_FREQ=ASNP_FREQ_CHS)]
  }

library(tictoc)
hap_origin=list()
for(chr in 1:22){
 tic(paste(pop,chr))
  ## read all snps
  ## read LD files
  ld <- fread(sprintf("%s/1KG/ld_1kG_%s_chr%s_popCellSNVs_nomono_hg38.ld.gz",DATA_DIR,pop,chr),h=T,stringsAsFactors = F, data.table = F)
  ld_b <- ld[,c("CHR_B","BP_B","SNP_B","CHR_A","BP_A","SNP_A","R2")] ## to make symmetric
  colnames(ld_b) <- colnames(ld)
  ld <- rbind(ld,ld_b)
  ld <- ld[ld$SNP_A%chin%archaic_IDs,]
  rm(ld_b)

  ## add intro info
  ld$posID <- ld$SNP_B
  ld <- merge(ld,aSNP_good[,.(posID,ORIGIN,ASNP_FREQ)],by="posID")
  ld$posID <- NULL
  ld$ORIGIN <- ifelse(is.na(ld$ORIGIN),"AMH",ld$ORIGIN)
  hap_origin[[chr]]=ld
  toc()
}
hap_origin=rbindlist(hap_origin)

hap_freq=hap_origin[,.(ASNP_FREQ=mean(ASNP_FREQ)),by=.(posID=SNP_A)]
asnps_test=merge(hap_freq,all_snps,by='posID')

#################################################################################################################
########   obtain global allele frequency, and frequency in population under consideration for all SNPs ########
#################################################################################################################

#############################################################################################
########   save the minimal set of data to reproduce the analysis and make the plots ########
#############################################################################################
file_name=sprintf('%s/tagaSNP_freqs/tag_aSNP_100kb_vs_eQTL_%s_%s_%s.tsv.gz',DAT_RESAMP_ASNP_DIR,pop,num_test,test_name)
asnps_table=asnps_test[(tag_SNP==TRUE | eQTL==TRUE) & DIST<1e5,.(posID, tag_SNP, eQTL, ASNP_FREQ, MEAN_MAF)]

fwrite(asnps_table,file=file_name,sep='\t')

stat_table=asnps_table[,.( POP=pop,
                            num=num_test,
                            eQTL_set=test_name,
                            DeltaFreq_eQTL=mean(ASNP_FREQ[eQTL])-mean(ASNP_FREQ[!eQTL]),
                            P_eQTL__wilcox=wilcox.test(ASNP_FREQ[eQTL],ASNP_FREQ[!eQTL])$p.value,
                            DeltaFreq_eQTL__adj_MAF=summary(lm(ASNP_FREQ~MEAN_MAF+eQTL))$coeff['eQTLTRUE',1],
                            P_eQTL__adj_MAF=anova(lm(ASNP_FREQ)~MEAN_MAF+eQTL))['eQTL',5],
                            P_eQTL_rank__adj_MAF=anova(lm(rank(ASNP_FREQ)~MEAN_MAF+eQTL))['eQTL',5],
                            )]

file_name_stats=sprintf('%s/users/Javier/results/tagaSNP_freqs/stats_tag_aSNP_100kb_vs_eQTL_%s_%s_%s.tsv',DAT_RESAMP_ASNP_DIR,pop,num_test,test_name)
fwrite(stat_table,file=file_name_stats,sep='\t')
