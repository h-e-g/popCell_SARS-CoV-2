
################################################################################
################################################################################
# File name: 4a2_aggregate_enrichment_results.R
# Author: J.MR, Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate enrichment results from 4a1 across all SNP lists
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


source(sprintf("%s/misc_plots.R",MISC_DIR))

snp_sets <- fread(sprintf("%s/All_eQTL_snpsSets.txt.gz",EQTL_DIR)) %>% select(set,snps)

allowed_celltypes=paste(c(lineage_order,celltype_order),collapse='|')
regex=sprintf('^(r?eQTL)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\2',set)]
snpSets[,state:=gsub(regex,'\\3',set)]
snpSets[,specificity:=gsub(regex,'\\4',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

####### resample PBS scores
resamples_PBS=dir(DAT_RESAMP_DIR,pattern='(YRI|CEU|CHS)_PBS_(.*)_NumTest(.*)_NumResamp10000.txt')

resamp_results=list()
for (i in resamples_PBS){
  cat(i,'\n')
  resamp_results[[i]]=fread(sprintf('%s/%s',DAT_RESAMP_DIR,i))
  resamp_results[[i]]=resamp_results[[i]][,.(RESAMP_NUM_SEL=mean(RESAMP_NUM_SEL),
      RESAMP_MEAN_SCORE=mean(RESAMP_MEAN_SCORE),
      FE_NUM_SEL=mean((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL)),
      lowerCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.025),
      upperCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.975),
      FE_MEAN_SCORE=mean(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,na.rm=T),
      lowerCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.025,na.rm=T),
      upperCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.975,na.rm=T))
      ,by=.(OBS_NUM_SEL,NUM_SEL_PVAL,OBS_MEAN_SCORE,MEAN_SCORE_PVAL,NUM_TARGET_EQTLS,NUM_RESAMP)]
}

resamp_results=rbindlist(resamp_results,idcol='snp_set_resamp')
REGEX=sprintf(".*/(CEU|YRI|CHS)_PBS_((r?eQTL_?(%s)?_*(NS|IAV|COV)?_?(shared|specific)?)_NumTest([0-9]+)_NumResamps10000.txt",allowed_celltypes)
resamp_results[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
resamp_results[,stat:='PBS']
resamp_results[,set:=gsub(REGEX,'\\3',snp_set_resamp)]
resamp_results[,num:=gsub(REGEX,'\\7',snp_set_resamp)]

resamp_results[,FDR_MEAN_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr'),by=stat]
resamp_results[,FDR_NUM_SEL:=p.adjust(NUM_SEL_PVAL,'fdr'),by=stat]

resamp_results=merge(resamp_results,unique(snpSets[,.(set, type, celltype, state, specificity)]),by='set')

#### exploration
# resamp_results[grepl('eQTL_(COV|IAV|NS)',set) & stat=='PBS',][,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL)][,FDR_NUM:=p.adjust(NUM_SEL_PVAL)][1:.N]
# resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',][,FDR:=p.adjust(MEAN_SCORE_PVAL,'fdr')][1:.N]

####### sets of eQTL tested (these are the ones that should be plotted in priority)
# consider all eQTL & reQTL at p<.01 in a given condition (sepArately for each condition),
# split reQTLs between shared and virus specific
resamp_results_USED=resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) & stat=='PBS',]
resamp_results_USED[,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr')]
resamp_results_USED[,FDR_NUM:=p.adjust(NUM_SEL_PVAL,'fdr')][1:.N]
fwrite(resamp_results_USED,file=sprintf("%s/result_resampling_allPOP_PBS_allSETs_10000Resamp.txt",DAT_RESAMP_DIR),sep='\t')

plot_resamp=function(resamp_results,stat=c('NUM','SCORE')){
  # a function that takes a subset of resampling results as input, and plots them.
  # eg. resamp_results[grepl('reQTL_(COV|IAV|NS)',set) & stat=='SCORE',]

  stat=match.arg(stat)
  df=resamp_results

  if(stat=='SCORE'){
    df_STAT=df[,.(point=FE_MEAN_SCORE,
        lowci=lowerCI_MEAN_SCORE,
        highci=upperCI_MEAN_SCORE,
        pval=MEAN_SCORE_PVAL),
        by=.(POP,type,celltype,state,specificity)]
      }
  if(stat=='NUM'){
    df_STAT=df[,.(point=FE_NUM_SEL,
        lowci=lowerCI_NUM_SEL,
        highci=upperCI_NUM_SEL,
        pval=NUM_SEL_PVAL),
        by=.(POP,type,celltype,state,specificity)]
        }
  df_STAT[,FDR:=p.adjust(pval,'fdr')]
  df_STAT[state=='NS',state:='b']
  df_STAT[,state2:=paste(celltype,state,specificity)]
  p <- ggplot(df_STAT,aes(x = state2,y=point,color=POP)) +
    geom_errorbar(aes(ymin=lowci, ymax=highci),width=0.2,position=position_dodge(width = 0.6)) +
    geom_point(aes(ymin=lowci, ymax=highci),position=position_dodge(width = 0.6)) +
    facet_grid(. ~ type, scales = "free") +
    scale_color_manual(values=c("#eba206","#71458d","#008000")) +
    geom_hline(aes(yintercept=1), colour="black",linetype="dotted") +
    xlab("") + ylab("Fold enrichment") + labs(color="Population") +
    theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
    p
}

FIG_DIR='../../../FIGURES/Fig4'

for (STAT in c('NUM','SCORE')){
    cat(STAT, 'PBS')
    # eQTL per condition
    pdf(sprintf('%s/Condition_eQTLs_PBS_%s.pdf',FIG_DIR,STAT),height=6.7*.4,width=7.2*.4)
    p <- plot_resamp(resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) ,],STAT)
    print(p)
    dev.off()

    ###### all eQTL
    # eQTL all celltypes
    pdf(sprintf('%s/Celltype_eQTLs_PBS_%s.pdf',FIG_DIR,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_(%s)$',celltype_order),set)  ,],STAT)
    print(p)
    dev.off()

    # reQTL all celltypes
    pdf(sprintf('%s/Celltype_reQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^reQTL_(%s)$',celltype_order),set) ,],STAT)
    print(p)
    dev.off()
}
