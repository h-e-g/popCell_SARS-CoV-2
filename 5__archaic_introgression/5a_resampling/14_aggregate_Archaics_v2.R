# definition of aSNPs and test for enrichments are in
#/evo_immuno_pop/users/Javier/scripts/
# definition archaic SNPs Table


options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
EIP=EVO_IMMUNO_POP_ZEUS

.libPaths(sprintf("%s/single_cell/resources/R_libs/4.1.0",EVO_IMMUNO_POP_ZEUS))

suppressMessages(library(data.table))
suppressMessages(library(tictoc))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
RES_DIR=sprintf("%s/single_cell/resources",EVO_IMMUNO_POP_ZEUS)

CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"

source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))
source(sprintf("%s/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",RES_DIR))


snpSets=fread(sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

# extract quantiles as we are reading
resampASNP_results_full=list()
resamples=dir(sprintf('%s/users/Javier/results/resampling_asnps_pruned/',EVO_IMMUNO_POP_ZEUS))
for (i in resamples){
  cat(i,'\n')
  resampASNP_results_full[[i]]=fread(sprintf('%s/users/Javier/results/resampling_asnps_pruned/%s',EVO_IMMUNO_POP_ZEUS,i))
  resampASNP_results_full[[i]]=resampASNP_results_full[[i]][,.(RESAMP_NUM_SEL=mean(RESAMP_NUM_SEL),
      RESAMP_MEAN_SCORE=mean(RESAMP_MEAN_SCORE),
      FE_NUM_SEL=(1+OBS_NUM_SEL)/mean(1+RESAMP_NUM_SEL),
      lowerCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.025),
      upperCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.975),
      FE_MEAN_SCORE=(1+OBS_NUM_SEL)/mean(1+RESAMP_NUM_SEL),
      lowerCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.025,na.rm=T),
      upperCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.975,na.rm=T))
      ,by=.(OBS_NUM_SEL,NUM_SEL_PVAL_ENRICHMENT,NUM_SEL_PVAL_DEPLETION,OBS_MEAN_SCORE,MEAN_SCORE_PVAL_ENRICHMENT,MEAN_SCORE_PVAL_DEPLETION,NUM_TARGET_EQTLS,NUM_RESAMP)]
}
resampASNP_results_full=rbindlist(resampASNP_results_full,idcol='snp_set_resamp')
REGEX=sprintf("(CEU|YRI|CHS)_(Any|Vindija33.19|Denisova|Archaic)_(FinalLenientaSNPs|FinalaSNPs|FinalAdaptiveaSNPs|FinalAdaptiveLenientaSNPs)_((r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?)_(.*)_NumResamps10000.txt",allowed_celltypes)
REGEX2="(CEU|YRI|CHS)_(Any|Vindija33.19|Denisova|Archaic)_(FinalLenientaSNPs|FinalaSNPs|FinalAdaptiveaSNPs|FinalAdaptiveLenientaSNPs)_(.*)_(adjlocalMAF.*)_NumResamps10000.txt"

resampASNP_results_full[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
resampASNP_results_full[,archaic:=gsub(REGEX,'\\2',snp_set_resamp)]
resampASNP_results_full[,aSNP:=gsub(REGEX,'\\3',snp_set_resamp)]
resampASNP_results_full[,set:=gsub(REGEX,'\\4',snp_set_resamp)]
resampASNP_results_full[,adj:=gsub(REGEX2,'\\5',snp_set_resamp)]
resampASNP_results_full=merge(resampASNP_results_full,unique(snpSets[,.(set, num, type, celltype, state, specificity)]),by='set')
resampASNP_results_full[specificity!='',type:=paste(type,'breakdown')]

fwrite(resampASNP_results_full,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/pruning/resampASNP_results_full.txt',EVO_IMMUNO_POP_ZEUS),sep='\t')




asnps=fread('/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt')
fwrite(asnps,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXa_asnps.tsv',EIP),sep='\t')


SuppTable_enrich_byCond=resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & state!='' & celltype=='' & type%in%c('eQTL','reQTL'),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]
SuppTable_enrich_byCond=SuppTable_enrich_byCond[,.(POP,type,state,aSNP_eQTL=OBS_NUM_SEL,FoldEnrich=FE_NUM_SEL,lowerCI_FE=lowerCI_NUM_SEL,upperCI_FE=upperCI_NUM_SEL,P=NUM_SEL_PVAL_ENRICHMENT,FDR)]
fwrite(SuppTable_enrich_byCond,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXb_enrich_byCond.tsv',EIP),sep='\t')

resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & type=='eQTL' & state!='',][order(lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]

resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & type=='eQTL' & state!='',][order(lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N][FDR<.1]
resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & type=='eQTL' & state!='' & celltype!='' & POP=='CEU',][order(lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N][FDR<.1]

resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & type=='eQTL' & state!='',][order(lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]

SuppTable_detail_enrich_byCelltype=resampASNP_results_full[aSNP=='FinalLenientaSNPs' & adj=='adjlocalMAF_LD_DIST_prunedBGD' & specificity=='' & archaic=='Any' & type=='eQTL' & state!='',][order(-lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]
SuppTable_detail_enrich_byCelltype=SuppTable_detail_enrich_byCelltype[,.(POP,type,celltype,state,aSNP_eQTL=OBS_NUM_SEL,FoldEnrich=FE_NUM_SEL,lowerCI_FE=lowerCI_NUM_SEL,upperCI_FE=upperCI_NUM_SEL,P=NUM_SEL_PVAL_ENRICHMENT,FDR)]
fwrite(SuppTable_detail_enrich_byCelltype,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXcdetail_enrich_byCelltype.tsv',EIP),sep='\t')


# per condition
stats[aSNP_set=='Any' & type=='eQTL' & celltype=='' & state!='' & specificity=='',][order(POP,state)]

SuppTable_FreqDiff_byCond=stats[aSNP_set=='Any' & type=='eQTL' & celltype=='' & state!='' & specificity=='',][order(POP,state)][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_FreqDiff_byCond=SuppTable_FreqDiff_byCond[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_FreqDiff_byCond,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXd_FreqDiff_byCond.tsv',EIP),sep='\t')

stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_detail_byCelltype=SuppTable_detail_byCelltype[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_detail_byCelltype,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXc_detailFreqDiff_byCelltype.tsv',EIP),sep='\t')
stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF_DIST_LD,'fdr')][FDR<.01,]


stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL__adj_MAF_DIST_LD,'fdr')][FDR<.01,]

stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')][FDR<.01,]


SuppTable_detail_byCelltype=stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_detail_byCelltype=SuppTable_detail_byCelltype[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_detail_byCelltype,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableSXXc_detailFreqDiff_byCelltype.tsv',EIP))
stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF_DIST_LD,'fdr')][FDR<.01,]


x=stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',]

x[,z:=abs(qnorm(P_eQTL__adj_MAF_DIST_LD/2))]
x[,se:=DeltaFreq_eQTL__adj_MAF_DIST_LD/z]
x[,delta_lower:=DeltaFreq_eQTL__adj_MAF_DIST_LD-2*se]
x[,delta_upper:=DeltaFreq_eQTL__adj_MAF_DIST_LD+2*se]
x[,FDR:=p.adjust(P_eQTL__adj_MAF_DIST_LD,'fdr')]
#x[FDR<.1 & P_eQTL_rank__adj_MAF_DIST_LD<.01,][order(-delta_lower),][1:10,]
x[ P_eQTL_rank__adj_MAF_DIST_LD<.01,][order(-delta_lower),][1:10,]


plot_resamp=function(resamp_results,stat=c('NUM','SCORE'),aSNP_set=c('FinalLenientaSNPs','FinalaSNPs','FinalAdaptiveaSNPs','FinalAdaptiveLenientaSNPs'),Archaic=c('Any','Vindija33.19','Denisova','Archaic'),ADJ){
  # a function that takes a subset of resampling results as input, and plots them.
  # eg. resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',]
  stat=match.arg(stat)
  aSNP_set=match.arg(aSNP_set)
  Archaic=match.arg(Archaic)
  df=resamp_results[aSNP==aSNP_set & archaic==Archaic & adj==ADJ,]
  if(stat=='SCORE'){
    df_STAT=df[,.(point=FE_MEAN_SCORE,
        lowci=lowerCI_MEAN_SCORE,
        highci=upperCI_MEAN_SCORE,
        pval=MEAN_SCORE_PVAL_ENRICHMENT),
        by=.(POP,type,celltype,state,specificity)]
      }
  if(stat=='NUM'){
    df_STAT=df[,.(point=FE_NUM_SEL,
        lowci=lowerCI_NUM_SEL,
        highci=upperCI_NUM_SEL,
        pval=NUM_SEL_PVAL_ENRICHMENT),
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


FIG_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/pruning',EVO_IMMUNO_POP_ZEUS)
for (aSNP_set in c('FinalLenientaSNPs','FinalAdaptiveLenientaSNPs')){
  for (ARCHAIC in c('Any','Vindija33.19','Archaic')){
    for (STAT in c('NUM')){
      for (ADJ in c('adjlocalMAF_DIST_prunedBGD','adjlocalMAF_LD_DIST_prunedBGD')){
        FIG_DIR_DETAIL=paste0(FIG_DIR,ifelse(grepl('LD',ADJ),ifelse(grepl('Adaptive',aSNP_set),'/Adaptive',''),'/noLD'))

        cat(STAT, ARCHAIC, aSNP_set,ADJ,'\n')
        # eQTL per condition
        pdf(sprintf('%s/Condition_eQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2*.4)
        p <- plot_resamp(resampASNP_results_full[grepl('eQTL_(COV|IAV|NS|shared)',set) & specificity!='stronger' ,],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        ###### all eQTL
        # eQTL all celltypes
        pdf(sprintf('%s/Celltype_eQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^eQTL_(%s)$',allowed_celltypes),set)  ,],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        # reQTL all celltypes

        pdf(sprintf('%s/Celltype_reQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^reQTL_(%s)$',allowed_celltypes),set) ,],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        ###### GW Signif eQTL

        pdf(sprintf('%s/Condition_GWsignif_eQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2*.4)
        p <- plot_resamp(resampASNP_results_full[grepl('GWsignif_(COV|IAV|NS)',set) ,],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        # eQTL all celltypes
        pdf(sprintf('%s/Celltype_GWsignif_eQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^eQTL_GWsignif_(%s)$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()
        # reQTL all celltypes
        pdf(sprintf('%s/Celltype_GWsignif_reQTLs_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^reQTL_GWsignif_(%s)$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        ###### all eQTL NS
        # eQTL NS all celltypes
        pdf(sprintf('%s/Celltype_eQTLs_NS_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^eQTL_(%s)__NS$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()


        ###### all eQTL COV
        # eQTL COV all celltypes
        pdf(sprintf('%s/Celltype_eQTLs_COV_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^eQTL_(%s)__COV$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        ###### all eQTL IAV
        # eQTL COV all celltypes
        pdf(sprintf('%s/Celltype_eQTLs_IAV_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^eQTL_(%s)__IAV$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()

        ###### all reQTL COV
        # eQTL COV all celltypes
        pdf(sprintf('%s/Celltype_reQTLs_COV_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^reQTL_(%s)__COV$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()


        ###### all reQTL IAV
        # reQTL IAV all celltypes
        pdf(sprintf('%s/Celltype_reQTLs_IAV_%s_%s_%s_%s.pdf',FIG_DIR_DETAIL,STAT,ARCHAIC,aSNP_set, ADJ),height=6.7*.4,width=7.2)
        p <- plot_resamp(resampASNP_results_full[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),set),],STAT, aSNP_set, ARCHAIC, ADJ)
        print(p)
        dev.off()
      }
    }
  }
}

################################################################
###### comparison of frequency between eQTL and non eQTL #######
################################################################

comparisons=dir(sprintf('%s/users/Javier/results/tagaSNP_freqs/stats',EVO_IMMUNO_POP_ZEUS))

stats=list()
for (i in comparisons){
  stats[[i]]=fread(sprintf('%s/users/Javier/results/tagaSNP_freqs/stats/%s',EVO_IMMUNO_POP_ZEUS,i))
}
stats=rbindlist(stats,idcol='comparison')
stats=merge(stats,unique(snpSets[,.(eQTL_set=set, num, type, celltype, state, specificity)]),by='eQTL_set')

stats[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),set),]

plot_FreqDiff=function(stats,archaic){
  stats0=stats[aSNP_set==archaic,]
  tags_SNPs=list()
  for (i in stats0[,comparison]){
    tags_SNPs[[i]]=fread(sprintf('%s/users/Javier/results/tagaSNP_freqs/%s.gz',EVO_IMMUNO_POP_ZEUS,gsub('stats_(.*)','\\1',i)))
  }
  tags_SNPs=rbindlist(tags_SNPs,idcol='comparison')
  tags_SNPs=merge(tags_SNPs,stats0,by='comparison')
  tags_SNPs[specificity!='',type:=paste(type,'breakdown')]
  tags_SNPs[state=='NS',state:='b']
  tags_SNPs[,state2:=paste(celltype,state,specificity)]
  tags_SNPs[,eQTL2:=factor(ifelse(eQTL,'eQTL','no eQTL'),c('no eQTL','eQTL'))]
  p <- ggplot(tags_SNPs,aes(x = state2,y=ASNP_FREQ,color=eQTL2,fill=POP)) +
    geom_violin(scale='width',size=.5) + geom_boxplot(fill='white',alpha=0.5,notch=TRUE,position=position_dodge(width=0.85),size=.3,outlier.size=.3) +
    facet_grid(POP ~ type, scales = "free") +
    scale_fill_manual(values=c("#eba206","#71458d","#008000")) +
    scale_color_manual(values=c('no eQTL'="#000000",'eQTL'="#e03342"))+
    xlab("") + ylab("archaic SNP frequency") + labs(color="Population") +
    theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + ylim(c(0,1))
    p

}

FIG_DIR_FREQ=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/FreqDiff',EVO_IMMUNO_POP_ZEUS)
dir.create(FIG_DIR_FREQ)
for (ARCHAIC in c('Any','Vindija33.19','Archaic')){
  cat(ARCHAIC,'\n')
  # eQTL per condition
  pdf(sprintf('%s/Condition_eQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*0.7,width=7.2*0.7)
  p <- plot_FreqDiff(stats[grepl('eQTL_(COV|IAV|NS|shared)',eQTL_set) & specificity!='stronger' ,],ARCHAIC)
  print(p)
  dev.off()
  stats[grepl('eQTL_(COV|IAV|NS|shared)',eQTL_set) & specificity!='stronger' & aSNP_set=='Any' ,.(aSNP_set,POP,type,celltype,DeltaFreq_eQTL,P_eQTL__wilcox,P_eQTL__adj_MAF_DIST_LD)][order(POP),]
  stats[grepl('eQTL_(COV|IAV|NS|shared)',eQTL_set) & specificity!='stronger' &  aSNP_set=='Vindija33.19',.(aSNP_set,POP,type,state,specificity,celltype,DeltaFreq_eQTL,P_eQTL__wilcox,P_eQTL__adj_MAF_DIST_LD)][order(POP),]

  ###### all eQTL
  # eQTL all celltypes
  pdf(sprintf('%s/Celltype_eQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^eQTL_(%s)$',allowed_celltypes),eQTL_set)  ,],ARCHAIC)
  print(p)
  dev.off()

  # stats[grepl(sprintf('^eQTL_(%s)$',allowed_celltypes),eQTL_set) & aSNP_set=='Any',.(aSNP_set,POP,type,celltype,DeltaFreq_eQTL,P_eQTL__wilcox,P_eQTL__adj_MAF_DIST_LD)][order(POP),]

  # reQTL all celltypes

  pdf(sprintf('%s/Celltype_reQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^reQTL_(%s)$',allowed_celltypes),eQTL_set) ,],ARCHAIC)
  print(p)
  dev.off()

  ###### GW Signif eQTL

  pdf(sprintf('%s/Condition_GWsignif_eQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*.7)
  p <- plot_FreqDiff(stats[grepl('GWsignif_(COV|IAV|NS)',eQTL_set) ,],ARCHAIC)
  print(p)
  dev.off()

  # eQTL all celltypes
  pdf(sprintf('%s/Celltype_GWsignif_eQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^eQTL_GWsignif_(%s)$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()
  # reQTL all celltypes
  pdf(sprintf('%s/Celltype_GWsignif_reQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^reQTL_GWsignif_(%s)$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()

  ###### all eQTL NS
  # eQTL NS all celltypes
  pdf(sprintf('%s/Celltype_eQTLs_NS_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^eQTL_(%s)__NS$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()


  ###### all eQTL COV
  # eQTL COV all celltypes
  pdf(sprintf('%s/Celltype_eQTLs_COV_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^eQTL_(%s)__COV$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()

  ###### all eQTL IAV
  # eQTL COV all celltypes
  pdf(sprintf('%s/Celltype_eQTLs_IAV_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^eQTL_(%s)__IAV$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()

  ###### all reQTL COV
  # eQTL COV all celltypes
  pdf(sprintf('%s/Celltype_reQTLs_COV_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^reQTL_(%s)__COV$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()


  ###### all reQTL IAV
  # reQTL IAV all celltypes
  pdf(sprintf('%s/Celltype_reQTLs_IAV_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*.7,width=7.2*1.8)
  p <- plot_FreqDiff(stats[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),eQTL_set),],ARCHAIC)
  print(p)
  dev.off()
  stats[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),eQTL_set) & aSNP_set=='Any',.(aSNP_set,POP,type,celltype,DeltaFreq_eQTL,P_eQTL__wilcox,P_eQTL__adj_MAF_DIST_LD)][order(POP),]
  stats[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),eQTL_set) & aSNP_set=='Vindija33.19',.(aSNP_set,POP,type,celltype,DeltaFreq_eQTL,P_eQTL__wilcox,P_eQTL__adj_MAF_DIST_LD)][order(POP),]


}
########################################
########### Final Figure 6a  ###########
########################################

# see /evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/pruning/Condition_eQTLs_NUM_Any_FinalLenientaSNPs_adjlocalMAF_LD_DIST_prunedBGD.pdf

########################################
########### Final Figure 6b  ###########
########################################

plot_FreqDiff=function(stats,archaic){
  stats0=stats[aSNP_set==archaic,]
  tags_SNPs=list()
  for (i in stats0[,comparison]){
    tags_SNPs[[i]]=fread(sprintf('%s/users/Javier/results/tagaSNP_freqs/%s.gz',EVO_IMMUNO_POP_ZEUS,gsub('stats_(.*)','\\1',i)))
  }
  tags_SNPs=rbindlist(tags_SNPs,idcol='comparison')
  tags_SNPs=merge(tags_SNPs,stats0,by='comparison')
  tags_SNPs[specificity!='',type:=paste(type,'breakdown')]
  tags_SNPs[state=='NS',state:='b']
  tags_SNPs[,state2:=paste(celltype,state,specificity)]
  tags_SNPs[,eQTL2:=factor(ifelse(eQTL,'eQTL','no eQTL'),c('no eQTL','eQTL'))]
  p <- ggplot(tags_SNPs,aes(x = state2,y=pmin(ASNP_FREQ,0.5),color=eQTL2,fill=POP)) +
    geom_violin(scale='width',size=.5) + geom_boxplot(fill='white',alpha=0.5,notch=TRUE,position=position_dodge(width=0.85),size=.3,outlier.size=.3) +
    facet_grid( ~ POP, scales = "free_x") +
    scale_fill_manual(values=c("#eba206","#71458d","#008000")) +
    scale_color_manual(values=c('no eQTL'="#000000",'eQTL'="#e03342"))+
    xlab("") + ylab("archaic SNP frequency") + labs(color="Population") +
    theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + ylim(c(0,.6))
    p

}

ARCHAIC='Any'
  pdf(sprintf('%s/02_Condition_eQTLs_%s_wide.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*0.4,width=7.2*0.5)
  p <- plot_FreqDiff(stats[grepl('eQTL_(COV|IAV|NS)',eQTL_set) & specificity=='' & type=='eQTL',],ARCHAIC)
  print(p)
  dev.off()

ARCHAIC='Any'
pdf(sprintf('%s/00_Condition_eQTLs_%s.pdf',FIG_DIR_FREQ,ARCHAIC),height=6.7*0.4,width=7.2*0.3)
p <- plot_FreqDiff(stats[grepl('eQTL_(COV|IAV|NS)',eQTL_set) & specificity=='' & type=='eQTL',],ARCHAIC)
print(p)
dev.off()
