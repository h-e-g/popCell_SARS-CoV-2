options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer))

EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
FIGURE_DIR = sprintf("%s/single_cell/project/pop_eQTL/figures/",EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
OUT_DIR = eQTL_DIR
SCRATCH="/pasteur/appa/scratch/mrotival/"

theme_set(theme_bw())
theme_update(
  text=element_text(family="serif",size=12),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))

RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5
#CHR=22
FDR_TH=0.01
CORRECT_CELLTYPE=FALSE
CONDITIONAL=FALSE
USE_JOINT=FALSE

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--run_name_celltype' | cmd[i]=='-r2' ){RUN_NAME_SECOND = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  # if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--correct_celltype' | cmd[i]=='-t' ){CORRECT_CELLTYPE = as.logical(cmd[i+1])}
  if (cmd[i]=='--conditional' | cmd[i]=='-c' ){CONDITIONAL = as.logical(cmd[i+1])}
  if (cmd[i]=='--add_joint_mapping' | cmd[i]=='-j' ){USE_JOINT = as.logical(cmd[i+1])}
}

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)
CONDITIONAL_char=ifelse(CONDITIONAL,'_zvalcond','')
CORRECT_CELLTYPE_char=ifelse(CORRECT_CELLTYPE,'_Nctadj','')
JOINT_char=ifelse(USE_JOINT,'_allcond_and_joint','_allcond')
################################################################################
####      merge results of eQTL mapping by lineage & cell types              ###
################################################################################

# RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
# RUN_NAME='lineage_condition_logFC__logFC__CellPropLineage_SVs_220409'
QTL_assoc_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
QTL_peak_lineage=fread(file=sprintf('%s/%s/dist_%s/independent_eQTLs%s_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
QTL_peak_lineage[,length(unique(snps))]
#[1] 9150
QTL_peak_lineage[,length(unique(gene))]
#[1] 5198


# RUN_NAME='celltype_condition___CellPropLineage_SVs_220409'
# RUN_NAME='celltype_condition_logFC__logFC__CellPropLineage_SVs_220409'
QTL_assoc_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',OUT_DIR,RUN_NAME_SECOND,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME_SECOND),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')


correct_factor=ifelse(CORRECT_CELLTYPE,length(cellstates),1)
bestSNP_list=unique(QTL_peak_lineage[,.(gene,snps,snp_score,round)])

QTL_assoc_remaining0=QTL_assoc_celltype[FDR<FDR_TH/correct_factor & converged==TRUE,]
QTL_assoc_remaining=QTL_assoc_remaining0

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
QTL_peak_celltype=merge(QTL_assoc_remaining0,bestSNP_list,by=c('gene','snps'))
QTL_peak_both=rbind(QTL_peak_lineage,QTL_peak_celltype)
QTL_peak_both[,length(unique(snps))]
QTL_peak_celltype[,length(unique(snps))]
QTL_peak_lineage[,length(unique(snps))]
QTL_peak_celltype[,length(unique(celltype)),keyby=.(snps,gene)][,.N,by=V1][,Pct:=round(N/sum(N),3)][1:.N]

QTL_peak_celltype[,celltype:=gsub('(.*)__(.*)','\\1',cellstate)]
QTL_peak_celltype[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
QTL_peak_both[,celltype:=gsub('(.*)__(.*)','\\1',cellstate)]
QTL_peak_both[,state:=gsub('(.*)__(.*)','\\2',cellstate)]

QTL_peak_celltype[,length(unique(snps)),by=gene][,.N,keyby=V1][,.(V1,N,N/sum(N))]
QTL_peak_celltype[,length(unique(state)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
QTL_peak_celltype[,length(unique(celltype)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]

fwrite(QTL_peak_celltype,file=sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltypeLevel_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
# write the lineage and celltype file in both folders
fwrite(QTL_peak_both,file=sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
fwrite(QTL_peak_both,file=sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME_SECOND,CIS_DIST_TEXT,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
