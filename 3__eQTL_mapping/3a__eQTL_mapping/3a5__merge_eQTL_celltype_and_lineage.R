################################################################################
################################################################################
# File name: 3a5__merge_eQTL_celltype_and_lineage.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate fine Mapped eQTLs from all lineage and celltype resolution
# same script is used for reQTLs
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))

# set defaults values
RUN_ID="RUN1"
RUN_NAME_LINEAGE="lineage_condition___CellPropLineage_SVs_RUN1"
RUN_NAME_CELLTYPE="celltype_condition___CellPropLineage_SVs_RUN1"
FDR_TH=0.01
# fixed
CIS_DIST=1e5

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_id' | cmd[i]=='-r' ){RUN_ID = cmd[i+1]} # ID of the run
	if (cmd[i]=='--run_name_lineage' | cmd[i]=='-r' ){RUN_NAME_LINEAGE = cmd[i+1]} # ID of the run with lineage stats
  if (cmd[i]=='--run_name_celltype' | cmd[i]=='-r2' ){RUN_NAME_CELLTYPE = cmd[i+1]} # ID of the run with celltypes stats
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
}

EQTL_DIR = "3_eQTL_mapping/SumStats"
SIGNIF_EQTL_DIR=sprintf("3_eQTL_mapping/SumStats/%s",RUN_ID)
dir.create(SIGNIF_EQTL_DIR)
COMBINED_EQTL_DIR_LINEAGE=sprintf('%s/%s/dist_%s',EQTL_DIR,RUN_NAME_LINEAGE,CIS_DIST_TEXT)
COMBINED_EQTL_DIR_CELLTYPE=sprintf('%s/%s/dist_%s',EQTL_DIR,RUN_NAME_CELLTYPE,CIS_DIST_TEXT)

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)

################################################################################
####      merge results of eQTL mapping by lineage & cell types              ###
################################################################################

QTL_assoc_lineage=fread(sprintf('%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',COMBINED_EQTL_DIR_LINEAGE))
QTL_peak_lineage=fread(file=sprintf('%s/%s/dist_%s/independent_eQTLs_allCond_%spctFDR.txt.gz',COMBINED_EQTL_DIR_LINEAGE,FDR_char))
### checks
#Number  of eQTL and genes
QTL_peak_lineage[,length(unique(snps))]
#[1] 9150
QTL_peak_lineage[,length(unique(gene))]
#[1] 5198

QTL_assoc_celltype=fread(sprintf('%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',COMBINED_EQTL_DIR_CELLTYPE))

bestSNP_list=unique(QTL_peak_lineage[,.(gene,snps,snp_score,round)])

# keep only FDR significant celtypes
QTL_assoc_celltype=QTL_assoc_celltype[FDR<FDR_TH & converged==TRUE,]
QTL_assoc_remaining=QTL_assoc_celltype

# progressively remove eQTLs that are not independant from QTL_assoc_celltype
cat(nrow(bestSNP_list),"independent eQTL identified - ",nrow(QTL_assoc_remaining),"SNP-gene association remaining\n")
explained_components=merge(QTL_assoc_remaining,bestSNP_list,by=c('gene','snps'))[,.(gene,snps,cellstate,top_component,FDR)]
QTL_assoc_remaining=QTL_assoc_remaining[!paste(gene,cellstate,top_component,sep='_')%chin%explained_components[,paste(gene,cellstate,top_component,sep='_')]]
cur_round=0
while(nrow(QTL_assoc_remaining)>0){
  cur_round=cur_round+1
  cat(nrow(bestSNP_list),"independent eQTL identified - ",nrow(QTL_assoc_remaining),"SNP-gene association remaining\n")
  SNP_score=QTL_assoc_remaining[,.(snp_score=sum(abs(pip_top_component))),by=.(gene, snps)]
  bestSNP=SNP_score[order(-snp_score),head(.SD,1),by=gene]
  bestSNP[,round:=cur_round]
  explained_components=merge(QTL_assoc_remaining,bestSNP,by=c('gene','snps'))[,.(gene,snps,cellstate,top_component,FDR)]
  QTL_assoc_remaining=QTL_assoc_remaining[!paste(gene,cellstate,top_component,sep='_')%chin%explained_components[,paste(gene,cellstate,top_component,sep='_')]]
  bestSNP_list=rbind(bestSNP_list,bestSNP)
}
QTL_peak_celltype=merge(QTL_assoc_celltype,bestSNP_list,by=c('gene','snps'))
QTL_peak_both=rbind(QTL_peak_lineage,QTL_peak_celltype)

### checks
#Number  of eQTL and genes at each level
QTL_peak_both[,length(unique(snps))]
QTL_peak_celltype[,length(unique(snps))]
QTL_peak_lineage[,length(unique(snps))]

QTL_peak_both[,celltype:=gsub('(.*)__(.*)','\\1',cellstate)]
QTL_peak_both[,state:=gsub('(.*)__(.*)','\\2',cellstate)]


tic('get gene_name for eQTL annotation')
feature_toUse=fread('3_eQTL_mapping/data/genes_to_use.tsv'),header=F)$V1

gtf <- rtracklayer::import(sprintf("references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
Feature_annot=as.data.table(gtf)[type=='gene',.(gene_id,gene_name,seqnames, start, end, strand,gene_type)]
Feature_annot[is.na(gene_name),gene_name:=gene_id]
Feature_annot[,seqnames:=as.character(seqnames)]
# for IAV_M and IAV_NS  we have two lines with the same genename (2 transcripts)
# we only keep the longest transcript
Feature_annot=Feature_annot[!(gene_id=="IAV_M" & end==784) & !(gene_id=="IAV_NS" & end==719),]
# we match the remaining features to feature_toUse
Feature_annot=Feature_annot[match(feature_toUse,gene_id),]
toc()

QTL_peak_both[,gene_name:=Feature_annot[match(gene,gene_id),gene_name]]

fwrite(QTL_peak_both,file=sprintf('%s/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR.txt.gz',SIGNIF_EQTL_DIR,FDR_char))
