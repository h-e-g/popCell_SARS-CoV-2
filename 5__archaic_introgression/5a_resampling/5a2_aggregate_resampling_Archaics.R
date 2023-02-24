
################################################################################
################################################################################
# File name: 5a2_aggregate_resampling_Archaics.R
# Author: J.MR, Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate enrichment results from 5a1 across all SNP lists
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

# declare usefule functions
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))

# load snp sets
snp_sets <- fread(sprintf("%s/All_eQTL_snpsSets.txt.gz",EQTL_DIR)) %>% select(set,snps)

allowed_celltypes=paste(c(lineage_order,celltype_order),collapse='|')
regex=sprintf('^(r?eQTL)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\2',set)]
snpSets[,state:=gsub(regex,'\\3',set)]
snpSets[,specificity:=gsub(regex,'\\4',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

####### read resulst from resampling aSNPs
resamples_aSNP=dir(DAT_RESAMP_ASNP_DIR,pattern='(YRI|CEU|CHS)_introgressed_(.*)_NumTest([0-9]+)_NumResamp10000.txt')

resamp_results=list()
for (i in resamples_aSNP){
  cat(i,'\n')
  resamp_results[[i]]=fread(sprintf('%s/%s',DAT_RESAMP_DIR,i))
  resamp_results[[i]]=resamp_results[[i]][,.(RESAMP_NUM_SEL=mean(RESAMP_NUM_SEL),
      FE_NUM_SEL=mean((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL)),
      lowerCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.025),
      upperCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.975),
      ,by=.(OBS_NUM_SEL,NUM_SEL_PVAL,NUM_SEL_PVAL_ENRICHMENT,NUM_SEL_PVAL_DEPLETION,NUM_TARGET_EQTLS,NUM_RESAMP)]
}

resamp_results=rbindlist(resamp_results,idcol='snp_set_resamp')
REGEX=sprintf(".*/(CEU|YRI|CHS)_introgressed_((r?eQTL_?(%s)?_*(NS|IAV|COV)?_?(shared|specific)?)_NumTest([0-9]+)_NumResamps10000.txt",allowed_celltypes)
resamp_results[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
resamp_results[,stat:='aSNP']
resamp_results[,set:=gsub(REGEX,'\\3',snp_set_resamp)]
resamp_results[,num:=gsub(REGEX,'\\7',snp_set_resamp)]

resamp_results[,FDR_NUM_SEL:=p.adjust(NUM_SEL_PVAL,'fdr'),by=stat]

resamp_results=merge(resamp_results,unique(snpSets[,.(set, type, celltype, state, specificity)]),by='set')

fwrite(resamp_results,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/pruning/resampASNP_results_full.txt',EVO_IMMUNO_POP_ZEUS),sep='\t')

SuppTable_enrich_byCond=resamp_results[specificity=='' & state!='' & celltype=='' & type%in%c('eQTL','reQTL'),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]
SuppTable_enrich_byCond=SuppTable_enrich_byCond[,.(POP,type,state,aSNP_eQTL=OBS_NUM_SEL,FoldEnrich=FE_NUM_SEL,lowerCI_FE=lowerCI_NUM_SEL,upperCI_FE=upperCI_NUM_SEL,P=NUM_SEL_PVAL_ENRICHMENT,FDR)]

fwrite(SuppTable_enrich_byCond,file=sprintf('%s/SupptableS8b_enrich_byCond.tsv',DAT_RESAMP_ASNP_DIR),sep='\t')

SuppTable_detail_enrich_byCelltype=resamp_results[type=='eQTL' & state!='',][order(-lowerCI_NUM_SEL),][,FDR:=p.adjust(NUM_SEL_PVAL_ENRICHMENT,'fdr')][1:.N]
SuppTable_detail_enrich_byCelltype=SuppTable_detail_enrich_byCelltype[,.(POP,type,celltype,state,aSNP_eQTL=OBS_NUM_SEL,FoldEnrich=FE_NUM_SEL,lowerCI_FE=lowerCI_NUM_SEL,upperCI_FE=upperCI_NUM_SEL,P=NUM_SEL_PVAL_ENRICHMENT,FDR)]
fwrite(SuppTable_detail_enrich_byCelltype,file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableS8c_detail_enrich_byCelltype.tsv',EIP),sep='\t')
