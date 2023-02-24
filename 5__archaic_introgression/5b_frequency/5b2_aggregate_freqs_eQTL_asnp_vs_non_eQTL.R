
################################################################################
################################################################################
# File name: 5b2_aggregate_freqs_eQTL_asnp_vs_non_eQTL.R
# Author: J.MR, Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate analysis of differences in frequency
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

######## aggregate 5b

comparisons=dir(sprintf('%s/tagaSNP_freqs/',DAT_RESAMP_ASNP_DIR),pattern='stats')

stats=list()
for (i in comparisons){
  stats[[i]]=fread(sprintf('%s/tagaSNP_freqs/%s',DAT_RESAMP_ASNP_DIR,i))
}
stats=rbindlist(stats,idcol='comparison')
stats=merge(stats,unique(snpSets[,.(eQTL_set=set, num, type, celltype, state, specificity)]),by='eQTL_set')

SuppTable_FreqDiff_byCond=stats[type=='eQTL' & celltype=='' & state!='' & specificity=='',][order(POP,state)][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_FreqDiff_byCond=SuppTable_FreqDiff_byCond[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_FreqDiff_byCond,file=sprintf('%s/TableS8/SupptableS8d_FreqDiff_byCond.tsv',FIG_DIR),sep='\t')

stats[aSNP_set=='Any' & type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_detail_byCelltype=SuppTable_detail_byCelltype[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_detail_byCelltype,file=sprintf('%s/TableS8/SupptableS8e_detailFreqDiff_byCelltype.tsv',FIG_DIR),sep='\t')

SuppTable_detail_byCelltype=stats[ type=='eQTL' & celltype!='' & state!='' & specificity=='',][,FDR:=p.adjust(P_eQTL_rank__adj_MAF,'fdr')]
SuppTable_detail_byCelltype=SuppTable_detail_byCelltype[1:.N,.(POP,type,celltype,state,DeltaFreq_eQTL__adj_MAF,P_eQTL_rank__adj_MAF,FDR)]
fwrite(SuppTable_detail_byCelltype,file=sprintf('%s/TableS8/SupptableS8e_detailFreqDiff_byCelltype.tsv',FIG_DIR))
