# test for enrichments of PBS nSL are in
#/evo_immuno_pop/users/Javier/scripts/0064_compute_enrichment_PBS_3factors.R
#/evo_immuno_pop/users/Javier/scripts/0064_compute_enrichment_PBS_3factors.sh
#/evo_immuno_pop/users/Javier/scripts/0031_compute_enrichment_nSL_3factors.R
#/evo_immuno_pop/users/Javier/scripts/0031_compute_enrichment_nSL_3factors.sh


options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'

.libPaths(sprintf("%s/single_cell/resources/R_libs/4.1.0",EVO_IMMUNO_POP_ZEUS))

suppressMessages(library(data.table))
suppressMessages(library(tictoc))
suppressMessages(library(rtracklayer))
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

####### resample PBS scores
resamples_PBS=paste0('resampling_newq99_30_06_2022/',dir(sprintf('%s/users/Javier/results/resampling_newq99_30_06_2022',EVO_IMMUNO_POP_ZEUS)))
####### resample nSL scores
resamples_nSL=paste0('resampling/',dir(sprintf('%s/users/Javier/results/resampling',EVO_IMMUNO_POP_ZEUS),pattern='nSL'))

resamp_results=list()
for (i in c(resamples_PBS,resamples_nSL)){
  cat(i,'\n')
  resamp_results[[i]]=fread(sprintf('%s/users/Javier/results/%s',EVO_IMMUNO_POP_ZEUS,i))
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
REGEX=sprintf(".*/(CEU|YRI|CHS)_(ABSnSL|PBS)_((r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?)(_NumTest[0-9]+)?_NumResamps10000.txt",allowed_celltypes)
resamp_results[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
resamp_results[,stat:=gsub(REGEX,'\\2',snp_set_resamp)]
resamp_results[,set:=gsub(REGEX,'\\3',snp_set_resamp)]
resamp_results[,num:=gsub('_NumTest','',gsub(REGEX,'\\9',snp_set_resamp))]

resamp_results[,FDR_MEAN_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr'),by=stat]
resamp_results[,FDR_NUM_SEL:=p.adjust(NUM_SEL_PVAL,'fdr'),by=stat]

resamp_results=merge(resamp_results,unique(snpSets[,.(set, type, celltype, state, specificity)]),by='set')

##### debug check what has not run
# all_tests=expand.grid(resamp_results[,unique(POP)],resamp_results[,unique(stat)],snpSets[,unique(set)])
# all_tests=data.table(all_tests)
# names(all_tests)=c('POP','stat','set')
# resamp_results_missing=merge(resamp_results[,.N,by=.(set,POP,stat)],all_tests,all.y=TRUE)[is.na(N),]
# resamp_results_missing=merge(resamp_results_missing,unique(snpSets[,.(set,num)]),by='set')
# fwrite(resamp_results_missing,file=sprintf('%s/users/Javier/toRerun_PBS.txt',EVO_IMMUNO_POP_ZEUS),sep='\t')


#### exploration
# resamp_results[grepl('eQTL_(COV|IAV|NS)',set) & stat=='PBS',][,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL)][,FDR_NUM:=p.adjust(NUM_SEL_PVAL)][1:.N]
# resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',][,FDR:=p.adjust(MEAN_SCORE_PVAL,'fdr')][1:.N]

####### TODO: sets of eQTL tested (these are the ones that should be plotted in priority)
# set 1=  consider all eQTL & reQTL at p<.01 in a given condition (seperately for each condition), split reQTLs between shared and virus specific PBS only
resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) & stat=='PBS' & specificity!='stronger',][,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr')][,FDR_NUM:=p.adjust(NUM_SEL_PVAL,'fdr')][1:.N]
# set 2=  consider GW singinficant eQTL & reQTL, PBS only (seperaltely for each condition)
resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',][,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr')][,FDR_NUM:=p.adjust(NUM_SEL_PVAL,'fdr')][1:.N]
# set 3=
resamp_results[grepl('(IAV)',set) & stat=='PBS' & specificity!='stronger' & POP=='CEU' & type=='eQTL' & celltype!='',][,FDR_SCORE:=p.adjust(MEAN_SCORE_PVAL)][,FDR_NUM:=p.adjust(NUM_SEL_PVAL)][1:.N]

# resamples=dir(sprintf('%s/users/Javier/results/resampling_newq99_30_06_2022',EVO_IMMUNO_POP_ZEUS))
#
# resamp_results_full=list()
# for (i in resamples){
#   cat(i,'\n')
#   resamp_results_full[[i]]=fread(sprintf('%s/users/Javier/results/resampling_newq99_30_06_2022/%s',EVO_IMMUNO_POP_ZEUS,i))
# }
# resamp_results_full=rbindlist(resamp_results_full,idcol='snp_set_resamp')
# REGEX=sprintf("(CEU|YRI|CHS)_(ABSnSL|PBS)_((r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?)(_NumTest[0-9]+)?_NumResamps10000.txt",allowed_celltypes)
# resamp_results_full[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
# resamp_results_full[,stat:=gsub(REGEX,'\\2',snp_set_resamp)]
# resamp_results_full[,set:=gsub(REGEX,'\\3',snp_set_resamp)]
# resamp_results_full[,num:=gsub('_NumTest','',gsub(REGEX,'\\9',snp_set_resamp))]
# resamp_results_full=merge(resamp_results_full,unique(snpSets[,.(set, type, celltype, state, specificity)]),by='set')


plot_resamp=function(resamp_results,stat=c('NUM','SCORE'),metric=c('PBS','ABSnSL')){
  # a function that takes a subset of resampling results as input, and plots them.
  # eg. resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',]

  stat=match.arg(stat)
  metric=match.arg(metric)
  df=resamp_results[stat==metric,]

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

FIG_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig5',EVO_IMMUNO_POP_ZEUS)

for (METRIC in c('PBS','ABSnSL')){
  for (STAT in c('NUM','SCORE')){
    cat(STAT, METRIC)
    # eQTL per condition
    pdf(sprintf('%s/Condition_eQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2*.4)
    p <- plot_resamp(resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) & specificity!='stronger' ,],STAT,METRIC)
    print(p)
    dev.off()

    ###### all eQTL
    # eQTL all celltypes
    pdf(sprintf('%s/Celltype_eQTLs_NUM_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_(%s)$',allowed_celltypes),set)  ,],STAT, METRIC)
    print(p)
    dev.off()

    # reQTL all celltypes

    pdf(sprintf('%s/Celltype_reQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^reQTL_(%s)$',allowed_celltypes),set) ,],STAT, METRIC)
    print(p)
    dev.off()

    ###### GW Signif eQTL

    pdf(sprintf('%s/Condition_GWsignif_eQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2*.4)
    p <- plot_resamp(resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) ,],STAT,METRIC)
    print(p)
    dev.off()

    # eQTL all celltypes
    pdf(sprintf('%s/Celltype_GWsignif_eQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_GWsignif_(%s)$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()
    # reQTL all celltypes
    pdf(sprintf('%s/Celltype_GWsignif_reQTLs_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^reQTL_GWsignif_(%s)$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()

    ###### all eQTL NS
    # eQTL NS all celltypes
    pdf(sprintf('%s/Celltype_eQTLs_NS_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_(%s)__NS$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()


    ###### all eQTL COV
    # eQTL COV all celltypes
    pdf(sprintf('%s/Celltype_eQTLs_COV_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_(%s)__COV$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()

    ###### all eQTL IAV
    # eQTL COV all celltypes
    pdf(sprintf('%s/Celltype_eQTLs_IAV_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^eQTL_(%s)__IAV$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()

    ###### all reQTL COV
    # eQTL COV all celltypes
    pdf(sprintf('%s/Celltype_reQTLs_COV_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^reQTL_(%s)__COV$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()


    ###### all reQTL IAV
    # reQTL IAV all celltypes

    pdf(sprintf('%s/Celltype_reQTLs_IAV_%s_%s.pdf',FIG_DIR,METRIC,STAT),height=6.7*.4,width=7.2)
    p <- plot_resamp(resamp_results[grepl(sprintf('^reQTL_(%s)__IAV$',allowed_celltypes),set),],STAT, METRIC)
    print(p)
    dev.off()
  }
}

#
# Je mets quelques idées ici sur la sélection.
#
# En Europe on a un enrichissement de PBS forts qui semble être le plus fort sur les eQTLs de la condition IAV et les eQTLs des Gamma Delta T cells. ça pourrait avoir du sens si les gamma delta étant donné que les  gamma delta semblent être une des premieres lignes de défense contre certains pathogènes.
# En particulier:
# ils ont un phenotype effecteur pré-programmé
# Les T gd résident dans les tissus et permettant une réponse immédiate à leur ligands.
# Il reconnaissent de nombreux pathogènes et sont parmi les premiers à secreter de l'IL17A pour initier la réponse immunitaire
# Ils semblent jouer un role dans la survie des nourissons exposé à la grippe. (= grosse pression de sélection)
#
# Cf cette review
# https://www.nature.com/articles/nri3384
# et ce papier.
# https://pubmed.ncbi.nlm.nih.gov/30170813/
