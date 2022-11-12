## load libs
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(stringr)
library(tidyr)
library(vroom)
library(tictoc)

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

# num_test <- 24
# num_resamples <- 20
# pop <- "CHS"

## snps sets
#snp_sets <- fread("data/snp_sets/snpSets_June2022.txt") %>% select(set,snps)
snp_sets <- fread("data/snp_sets/snpSets_topSNPs_nov2022.txt.gz") %>% select(set,snps)
# length(unique(snp_sets$set))
# [1] 256
test_name <- unique(snp_sets$set)[num_test]

cat("PBS", pop, 'test nb:',num_test,':',test_name,'\n')

## original mapping file
tot_snps <- fread("results/snps_ids/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_ID_SNPS_REF_ALT.txt") %>% select(ID,snps)
# colnames(tot_snps)[3:4] <- c("REF_ORI","ALT_ORI")
gc()

## read files with factors to account for
all_snps <- fread(paste0("results/bins/",pop,"_ALLCHR_popCellSNVs_nomono_maf0.05_hg38_3BINS_PBS.txt.gz")) %>% select(ID,MEAN_MAF_BIN_0.01,DIST_BIN,LDscore_BIN)

## grab only snps from original mapping file
all_snps <- all_snps[all_snps$ID %in% tot_snps$ID,]

## make bin
all_snps <- all_snps %>% mutate(BIN=paste0(LDscore_BIN,"_",MEAN_MAF_BIN_0.01,"_",DIST_BIN)) %>% select(ID,BIN)

## get eqtls ID
snp_sets <- snp_sets[snp_sets$set==test_name,]
snp_sets <- left_join(snp_sets,tot_snps,by="snps")
eqtls <- unique(snp_sets$ID)
rm(tot_snps,snp_sets)
gc()

## read PBS
sel_stat <- fread(paste0("results/pbs/",pop,"_ALLCHR_popCellSNVs_nomono_hg38.pbs.txt.gz")) %>% select(ID,PBS)
quant99 <- quantile(sel_stat$PBS,0.99,na.rm=T)
gc()

## add PBS scores to all_snps
all_snps <- left_join(all_snps,sel_stat,by="ID") %>% na.omit()
all_snps$SELECTED <- ifelse(all_snps$PBS>=quant99,T,F)
rm(sel_stat);gc()

## check if eQTLs are present in file
eqtls_final <- eqtls[(eqtls %in% all_snps$ID)]
rm(eqtls)
gc()

## change columns to SCORE and SELECTED
colnames(all_snps) <- c("ID","BIN","SCORE","SELECTED")

## load function
source("scripts/resample_SNPs_matchedbyBIN_SCORE_SELECTED_function.R")

## resample
res <- resample_ids(target_snps = eqtls_final, all_snps_bins = all_snps, num_resamples = num_resamples, rm_target_from_all_snps = T, use_data.table=TRUE)
res <- as_tibble(res)

## save
write.table(res,
            file=paste0("results/resampling_newq99_30_06_2022/",pop,"_PBS_",test_name,"_NumTest",num_test,"_NumResamps",format(num_resamples,scientific=F),".txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")
