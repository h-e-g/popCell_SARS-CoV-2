## resample a vector of SNPs based on a DF with BINS

## load libs
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

## define function
resample_ids <- function(target_snps=NULL,all_snps_bins=NULL,rm_target_from_all_snps=F,num_resamples=1000,use_data.table=FALSE){

  # target_snps <- eqtls_final
  # all_snps_bins <- all_snps

  ## target_snps -- vector with unique elements (i.e. SNPs to test)
  ## all_snps_bins -- data frame with all SNPs (i.e background SNPs)
  ## columns defined as 1) ID (same as target_snps)
  ##                    2) BIN (factors that will be used to match and resample)
  ##                    3) SCORE (selection score)
  ##                    4) SELECTED (SNP is selected/non-selected -- boolean of T/F)
  ## num_resamples -- number of resamples
  ## rm_target_from_all_snps -- remove the target SNPs from the background SNP list (i.e. you cannot resample those)
  ## use_data.table -- reimplementation of SNP resampling based on datatable. provides significant speed up.

  ## get observed value
  obs <- left_join(data.frame(ID=target_snps),all_snps_bins,by="ID")
  obs_sel <- sum(obs$SELECTED)
  obs_sel_score <- mean(obs$SCORE)
  obs_bin <- sort(obs$BIN)
  freq_bin <- as_tibble(as.data.frame(table(obs_bin)))
  colnames(freq_bin)[1] <- "BIN"
  rm(obs)

  ## get only relevant bins -- to make all_snps_bins smaller
  all_snps_bins <- all_snps_bins[all_snps_bins$BIN %in% unique(obs_bin),]

  ## remove target snps from background snps
  if (rm_target_from_all_snps==T){
    all_snps_bins <- all_snps_bins[!(all_snps_bins$ID %in% target_snps),]
  }

  ## resample
  resamp_num_sel_vec <- c()
  resamp_sel_score_vec <- c()

  ## add freq_bin info to data.table
  if(use_data.table){
    all_snps_bins=merge(all_snps_bins,freq_bin,by='BIN')
  }

  for (i in 1:num_resamples){

    ## print
    print(paste(i,"/",num_resamples))

    ## resample SNPs matched by bin without replacement
if(use_data.table){
    resamp_df <- all_snps_bins[,.SD[sample(1:.N,unique(Freq),replace=F),.(SCORE,SELECTED)],by=BIN]
  }else{
    resamp_df <- all_snps_bins %>%
      nest(-BIN) %>%
      left_join(freq_bin, by = "BIN") %>%
      mutate(Sample = map2(data, Freq, sample_n)) %>%
      unnest(Sample) %>%
      select(SCORE,SELECTED)
    }

    resamp_num_sel_vec <- c(resamp_num_sel_vec,sum(resamp_df$SELECTED))
    resamp_sel_score_vec <- c(resamp_sel_score_vec,mean(resamp_df$SCORE))

  }

  ## compute one-sided P-values
  pval_num_sel <- sum(obs_sel <= resamp_num_sel_vec)/num_resamples
  pval_mean_sel <- sum(obs_sel_score <= resamp_sel_score_vec)/num_resamples

  pval_num_sel_dep <- sum(obs_sel >= resamp_num_sel_vec)/num_resamples
  pval_mean_sel_dep <- sum(obs_sel_score >= resamp_sel_score_vec)/num_resamples


  ## make data frame to output
  res <- data.frame(RESAMP_NUM_SEL=resamp_num_sel_vec,
                    RESAMP_MEAN_SCORE=resamp_sel_score_vec,

                    OBS_NUM_SEL=obs_sel,
                    OBS_MEAN_SCORE=obs_sel_score,

                    NUM_SEL_PVAL_ENRICHMENT=pval_num_sel,
                    MEAN_SCORE_PVAL_ENRICHMENT=pval_mean_sel,


                    NUM_SEL_PVAL_DEPLETION=pval_num_sel_dep,
                    MEAN_SCORE_PVAL_DEPLETION=pval_mean_sel_dep,

                    NUM_RESAMP=num_resamples,
                    NUM_TARGET_EQTLS=length(target_snps))
  return(res)
  # hist(res$RESAMP_MEAN_SCORE,main=paste0("P=",res$MEAN_SCORE_PVAL[1]));abline(v=res$OBS_MEAN_SCORE[1],col="red")
  # hist(res$RESAMP_NUM_SEL,main=paste0("P=",res$NUM_SEL_PVAL[1]));abline(v=res$OBS_NUM_SEL[1],col="red")

}
