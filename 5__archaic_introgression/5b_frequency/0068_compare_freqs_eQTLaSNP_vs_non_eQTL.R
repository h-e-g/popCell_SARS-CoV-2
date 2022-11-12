EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
setwd(sprintf('%s/users/Javier/',EIP))

library(data.table)
library(dplyr)

library(ggplot2)
RES_DIR=sprintf("%s/single_cell/resources",EIP)
source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))

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
pop=as.character(temp[2])
archaic=as.character(temp[3])

# num_test <- 109
# pop <- "CEU"
# archaic <- "Archaic"

## snps sets
snp_sets <- fread("data/snp_sets/snpSets_June2022.txt") %>% select(set,snps)
test_name <- unique(snp_sets$set)[num_test]

# cat("Archaic Lenient", archaic, ',', pop, 'test nb:',num_test,':',test_name,"- adjusted on", MAF,'MAF,',ifelse(LD,'LDscore,',''), ifelse(DIST,' and DIST',''), ifelse(PRUNE_BGD,' pruning background',''),'\n')


## read files with factors to account for in this case all SNPs
all_snps <- fread(paste0("results/bins/",pop,"_ALLCHR_popCellSNVs_nomono_maf0.05_hg38_3BINS_4lenientaSNPs_withGlobalMAF.txt.gz"))
all_snps <- all_snps %>% select(ID,MAF_BIN_0.01,DIST_BIN,LDscore_BIN,MEAN_MAF_BIN_0.01)
## define the univers of SNPs to consider
tot_snps <- fread("results/snps_ids/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_ID_SNPS_REF_ALT.txt") %>% select(ID,snps)
all_snps = left_join(tot_snps,all_snps)


## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,]
snp_sets <- left_join(snp_sets,tot_snps,by="snps")
rm(tot_snps)
gc()

eqtls <- unique(snp_sets$ID)
eqtls_snps <- unique(snp_sets$snps)

# PRUNE SET OF EQTLs
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


LD_pruned=fread(sprintf("%s/users/Javier/VCF/prune/tagSNPs_%s_maf0.05pct_LD0.8_allCHR.prune.in",EIP,pop),header=F)$V1
all_snps[,tag_SNP:=ID%chin%LD_pruned]
all_snps[,eQTL:=ID%in%eqtls]


##############################################
########   annotate introgressed SNPs ########
##############################################


if (archaic=="Vindija33.19"){
  archaic2 <- "NEAND"
  } else if (archaic=="Denisova") {
  archaic2 <- "DENI"
} else if (archaic=="Archaic" | archaic=='Any') {
  archaic2 <- toupper(archaic)
}
asnps <- fread(paste0("results/asnps/Final_lenient_Adaptive_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt"))
asnps <- asnps[Population==pop & HAP_ORIGIN==archaic2,]
archaic_IDs=asnps[,ID]



all_snps = left_join(all_snps,asnps[,.(ID,HAP_ORIGIN)])

########################################################################################
########   obtain frequency of introgressed allele the pop under consideration ########
########################################################################################


aSNP_good=fread('/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt')

library(tictoc)
hap_origin=list()
for(chr in 1:22){
 tic(paste(pop,chr))
  ## read all snps
  ## read LD files
  ld <- fread(paste0(EIP,"/users/Javier/results/ld/",pop,"_chr",chr,"_popCellSNVs_nomono_hg38.ld.gz"),h=T,stringsAsFactors = F, data.table = F)
  ld_b <- ld[,c("CHR_B","BP_B","SNP_B","CHR_A","BP_A","SNP_A","R2")] ## to make symmetric
  colnames(ld_b) <- colnames(ld)
  ld <- rbind(ld,ld_b)
  ld <- ld[ld$SNP_A%chin%archaic_IDs,]
  rm(ld_b)

  ## add intro info
  ld$ID <- ld$SNP_B
  ld <- merge(ld,aSNP_good[Population==pop,c("ID","ORIGIN","ASNP_FREQ")],by="ID")
  ld$ID <- NULL
  ld$ORIGIN <- ifelse(is.na(ld$ORIGIN),"AMH",ld$ORIGIN)
  hap_origin[[chr]]=ld
  toc()
}
hap_origin=rbindlist(hap_origin)

hap_freq=hap_origin[,.(ASNP_FREQ=mean(ASNP_FREQ)),by=.(ID=SNP_A)]
asnps_test=merge(hap_freq,all_snps,by='ID')


#################################################################################################################
########   obtain global allele frequency, and frequency in populationb under consideration for all SNPs ########
#################################################################################################################

asnps_test[,MEAN_MAF:=as.numeric(gsub('.(0.[0-9]+),(0.[0-9]+).','\\1',as.character(MEAN_MAF_BIN_0.01)))+0.005]
asnps_test[,MEAN_MAF_BIN0.05:=cut(MEAN_MAF,seq(0,0.5,by=0.05))]
asnps_test[,MAF:=as.numeric(gsub('.(0.[0-9]+),(0.[0-9]+).','\\1',as.character(MAF_BIN_0.01)))+0.005]


#############################################################################################
########   save the minimal set of data to reproduce the analysis and make the plots ########
#############################################################################################
file_name=sprintf('%s/users/Javier/results/tagaSNP_freqs/tag_aSNP_100kb_vs_eQTL__%s_%s_%s_%s.tsv.gz',EIP,pop,archaic,num_test,test_name)
asnps_table=asnps_test[(tag_SNP==TRUE | eQTL==TRUE) & DIST_BIN!='(1e+05,Inf]',.(ID, tag_SNP, eQTL, ASNP_FREQ, MEAN_MAF, MAF, DIST_BIN, LDscore_BIN)]

fwrite(asnps_table,file=file_name,sep='\t')

stat_table=asnps_table[,.( POP=pop,
                            aSNP_set=archaic,
                            num=num_test,
                            eQTL_set=test_name,
                            DeltaFreq_eQTL=mean(ASNP_FREQ[eQTL])-mean(ASNP_FREQ[!eQTL]),
                            P_eQTL__wilcox=wilcox.test(ASNP_FREQ[eQTL],ASNP_FREQ[!eQTL])$p.value,
                            DeltaFreq_eQTL__adj_MAF=summary(lm(ASNP_FREQ~MEAN_MAF+eQTL))$coeff['eQTLTRUE',1],
                            P_eQTL__adj_MAF=anova(lm(rank(ASNP_FREQ)~MEAN_MAF+eQTL))['eQTL',5],
                            P_eQTL_rank__adj_MAF=anova(lm(ASNP_FREQ~MEAN_MAF+eQTL))['eQTL',5],
                            DeltaFreq_eQTL__adj_MAF_DIST_LD=summary(lm(ASNP_FREQ~MEAN_MAF++DIST_BIN+LDscore_BIN+eQTL))$coeff['eQTLTRUE',1],
                            P_eQTL__adj_MAF_DIST_LD=anova(lm(ASNP_FREQ~MEAN_MAF+DIST_BIN+LDscore_BIN+eQTL))['eQTL',5],
                            P_eQTL_rank__adj_MAF_DIST_LD=anova(lm(rank(ASNP_FREQ)~MEAN_MAF+DIST_BIN+LDscore_BIN+eQTL))['eQTL',5]
                            )]

file_name_stats=sprintf('%s/users/Javier/results/tagaSNP_freqs/stats/stats_tag_aSNP_100kb_vs_eQTL__%s_%s_%s_%s.tsv',EIP,pop,archaic,num_test,test_name)
fwrite(stat_table,file=file_name_stats,sep='\t')

file_name_plots=sprintf('%s/users/Javier/results/tagaSNP_freqs/plots/plots_tag_aSNP_100kb_vs_eQTL__%s_%s_%s_%s.pdf',EIP,pop,archaic,num_test,test_name)
pdf(file_name_plots,width=2.2,height=4)
p <- ggplot(asnps_table,aes(x=eQTL,y=ASNP_FREQ))+geom_violin(scale='width',fill=setNames(color_populations,c('YRI','CEU','CHS'))[pop])+geom_boxplot(notch=TRUE,fill='white',alpha=0.5)+theme_yann()
print(p)
dev.off()


# asnps_test[(tag_SNP==TRUE | eQTL==TRUE),.(.N,Freq_eQTL_aSNP=mean(ASNP_FREQ[eQTL==TRUE]),
#                                                 Freq_neQTL_aSNP=mean(ASNP_FREQ[eQTL==FALSE]),
#                                                 P_wilcox=wilcox.test(ASNP_FREQ[eQTL==TRUE],ASNP_FREQ[eQTL==FALSE])$p.value)]
#
# asnps_test[(tag_SNP==TRUE | eQTL==TRUE),.(.N,Freq_eQTL_aSNP=mean(ASNP_FREQ[eQTL==TRUE]),
#                                                 Freq_neQTL_aSNP=mean(ASNP_FREQ[eQTL==FALSE]),
#                                                 ifelse(sum(eQTL==FALSE)>0 & sum(eQTL==TRUE)>0 ,wilcox.test(ASNP_FREQ[eQTL==TRUE],ASNP_FREQ[eQTL==FALSE])$p.value,NA)),by=MEAN_MAF_BIN0.05]
#
# res_wilcox=asnps_test[(tag_SNP==TRUE | eQTL==TRUE),.(.N, MEAN_MAF=mean(MEAN_MAF),
#                                                 Freq_eQTL_aSNP=mean(ASNP_FREQ[eQTL==TRUE]),
#                                                 Freq_neQTL_aSNP=mean(ASNP_FREQ[eQTL==FALSE]),
#                                                 P_wilcox=ifelse(sum(eQTL==FALSE)>0 & sum(eQTL==TRUE)>0 ,wilcox.test(ASNP_FREQ[eQTL==TRUE],ASNP_FREQ[eQTL==FALSE])$p.value,-1)),by=MEAN_MAF_BIN0.05]
#
# asnps_test[(tag_SNP==TRUE | eQTL==TRUE) & DIST_BIN!='(1e+05,Inf]',anova(lm(ASNP_FREQ~MEAN_MAF_BIN0.05+eQTL))]
# asnps_test[(tag_SNP==TRUE | eQTL==TRUE) & DIST_BIN!='(1e+05,Inf]',anova(lm(rank(ASNP_FREQ)~rank(MEAN_MAF)+eQTL+DIST_BIN+LDscore_BIN))]
# asnps_test[(tag_SNP==TRUE | eQTL==TRUE) & DIST_BIN!='(1e+05,Inf]',anova(lm(rank(ASNP_FREQ)~rank(MEAN_MAF)+eQTL+DIST_BIN+LDscore_BIN))]


# asnps_test[(tag_SNP==TRUE),summary(lm(eQTL~MEAN_MAF))]
# asnps_test[tag_SNP==TRUE,summary(lm(eQTL~MEAN_MAF+MAF))]
# asnps_test[tag_SNP==TRUE,summary(glm(eQTL~I(MEAN_MAF*100)+I(MAF*100)),family=binomial)]
#
# asnps_test[tag_SNP==TRUE,summary(glm(eQTL~MEAN_MAF+MAF),family=binomial)]
