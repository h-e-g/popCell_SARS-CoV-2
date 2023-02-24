#!/bin/bash
################################################################################
################################################################################
# File name: 6a_COVIDloci_enrichments.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: test enrichment of eQTL and reQTLS in COVID risk loci
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


# default parameters
NSAMPLE=10000
SET_NUM=1

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--set' | cmd[i]=='-s' ){SET_NUM = cmd[i+1]} # ID of the set to test
  if (cmd[i]=='--nsamp' | cmd[i]=='-n' ){NSAMPLE = cmd[i+1]} # ID of the set to test
}

###############################################
########### SNP informations genome wide ######
###############################################
tic('loding SNP info')
SNP_info=getMap(annotate=TRUE)
SNP_info[,AF_global:=(ALT_freq_AFB+ALT_freq_ASH+ALT_freq_CEU)/3]
toc()

###############################################
###############################################
########### resampling analyses ###############
###############################################
###############################################

########### resampling functions
#### step 1: resample SNPs and extract relevant informations

resample_matched_on_MAF_DIST=function(snpSet,NSAMP=1000){
  SNP_metadata=SNP_info[max_MAF>.05 & !is.na(covid_A2_pval*covid_B2_pval*covid_C2_pval) & DistGene<1e5,]
  SNP_metadata[,is_eQTL:=ID%chin%snpSet]
  SNP_metadata[,BIN:=paste(cut(pmin(AF_global,1-AF_global),seq(0,1,by=.01)),cut(DistGene/1000,c(-1,0,1,5,10,20,50,100,Inf)))]
  resamp_eQTL=list()
  resamp_eQTL[[1]]=SNP_metadata[is_eQTL==TRUE,.(BIN,covid_A2_pval,covid_B2_pval,covid_C2_pval)]

    for (i in 1:NSAMP){
      ## print
      tic(paste(i,"/",NSAMP))
      resamp_eQTL[[i+1]] <- SNP_metadata[,.SD[sample(1:.N,sum(is_eQTL),replace=F),.(covid_A2_pval,covid_B2_pval,covid_C2_pval)],by=BIN]
      toc()
    }
    resamp_eQTL=rbindlist(resamp_eQTL,idcol='resamp')
    resamp_eQTL[,resamp:=resamp-1]
  resamp_eQTL
}

#### step 2: compute Fold enrichments and Pvalues from a set of resampled SNPs.
generate_result_table=function(resamp_eQTL){
  # we add 1 pseudo count to the number of significant pvalues to avoid numerical issues
  resamp_eQTL_PctSignif=as.data.table(resamp_eQTL)[,.(
              Pct_A2_5pct=mean(c(1,covid_A2_pval<0.05),na.rm=T),
              Pct_A2_1pct=mean(c(1,covid_A2_pval<0.01),na.rm=T),
              Pct_A2_01pct=mean(c(1,covid_A2_pval<0.001),na.rm=T),
              Pct_A2_001pct=mean(c(1,covid_A2_pval<0.0001),na.rm=T),
              Pct_B2_5pct=mean(c(1,covid_B2_pval<0.05),na.rm=T),
              Pct_B2_1pct=mean(c(1,covid_B2_pval<0.01),na.rm=T),
              Pct_B2_01pct=mean(c(1,covid_B2_pval<0.001),na.rm=T),
              Pct_B2_001pct=mean(c(1,covid_B2_pval<0.0001),na.rm=T),
              Pct_C2_5pct=mean(c(1,covid_C2_pval<0.05),na.rm=T),
              Pct_C2_1pct=mean(c(1,covid_C2_pval<0.01),na.rm=T),
              Pct_C2_01pct=mean(c(1,covid_C2_pval<0.001),na.rm=T),
              Pct_C2_001pct=mean(c(1,covid_C2_pval<0.0001),na.rm=T)),by=resamp]


  resamp_eQTL_PctSignif=melt(resamp_eQTL_PctSignif,id.vars='resamp')
  resamp_eQTL_PctSignif[,FE:=value[resamp==0]/value,by=variable]
  resamp_eQTL_PctSignif[,threshold:=gsub('Pct_([A-C]2)_(.*)pct','\\2',variable)]
  resamp_eQTL_PctSignif[,threshold:=factor(paste0('0.0',threshold),c('0.05','0.01','0.001','0.0001'))]
  resamp_eQTL_PctSignif[,GWAS_type:=factor(gsub('Pct_([A-C]2)_(.*)pct','\\1',variable),c('C2','B2','A2'))]
  resamp_eQTL_PctSignif[,GWAS_type_full:=case_when(GWAS_type=='C2'~'reported',GWAS_type=='B2'~'hospitalized',GWAS_type=='A2'~'critical')]
  resamp_eQTL_PctSignif[,GWAS_type_full:=factor(GWAS_type_full,c('reported','hospitalized','critical'))]

result_table=resamp_eQTL_PctSignif[resamp!=0,.(FE_q025=quantile(FE,.025,na.rm=T),
                                              FE=mean(FE,na.rm=T),
                                              FE_q975=quantile(FE,.975,na.rm=T),
                                              pvalue=pnorm(0, mean(log(FE),na.rm=T),
                                                  sd(log(FE),na.rm=T),lower=T)/2),by=.(threshold,GWAS_type_full,GWAS_type)]

  for( i in 1:result_table[,.N]){
    colname=paste0('Pct_',result_table[i,GWAS_type],'_',gsub('0.0','',result_table[i,threshold]),'pct')
    observed=resamp_eQTL_PctSignif[resamp==0,get(colname)]
    result_table[i,P_resamp:=resamp_eQTL_PctSignif[,mean(get(colname)>=observed)]]
  }
  result_table[1:.N]
}

#### step 3: make an enrichment plot from an OR table
plot_enrichment=function(Enrich_table,DIR,pname){
  Estim_PctSignif=Enrich_table[is.finite(FE),][,.(FE=mean(FE),lower=quantile(FE,.025),upper=quantile(FE,.975)),by=.(threshold,GWAS_type_full)]
  color_COVID=c(critical="#DD5555", hospitalized="#E68686", reported="#F1BABA")

  pdf(sprintf("%s/%s",DIR,pname),width=4,height=6)
  p <- ggplot(Enrich_table)+geom_point(aes(x=threshold,y=FE,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5))
  p <- p + geom_errorbar(aes(x=threshold,y=FE,ymin=FE_q025,ymax=FE_q975,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5)) + theme_yann()
  p <- p + geom_hline(yintercept=1,col='grey') + xlab('p-value threshold')+ ylab('Fold Enrichment of eQTLs in significant associations') + ylim(c(0,15))
  p <- p + scale_color_manual(values=color_COVID)
  print(p)
  dev.off()
}

snpSET=snpSets[num==SET_NUM,snps]
set.NAME=snpSets[num==SET_NUM,make.names(unique(set))]

# eQTL snps
snpSET=snpSET[!is.na(SNP_info[match(snpSET,ID),covid_A2_pval*covid_B2_pval*covid_C2_pval])]

resamp_eQTL=resample_matched_on_MAF_DIST(snpSET,NSAMP=NSAMPLE)

fwrite(resamp_eQTL,file=sprintf("%s/resamp_%sx_%s.tsv.gz",DAT_RESAMP_GWAS_DIR,NSAMPLE,set.NAME),sep='\t')
resamp_eQTL_Enrichment=generate_result_table(resamp_eQTL)

fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/resamp_COVID19_Enrichment_%sx_%s.tsv.gz",DAT_RESAMP_GWAS_DIR,NSAMPLE,set.NAME),sep='\t')

plot_enrichment(resamp_eQTL_Enrichment,DAT_RESAMP_GWAS_DIR,sprintf('resamp_COVID19_Enrichment_%sx_%s_plot.pdf',NSAMPLE,set.NAME))

q('no')
