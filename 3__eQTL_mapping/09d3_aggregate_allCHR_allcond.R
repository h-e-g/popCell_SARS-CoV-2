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

RUN_NAME="/lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5
#CHR=22
FDR_TH=0.01
CORRECT_CELLTYPE=FALSE
CONDITIONAL=FALSE

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  # if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--correct_celltype' | cmd[i]=='-t' ){CORRECT_CELLTYPE = as.logical(cmd[i+1])}
  if (cmd[i]=='--conditional' | cmd[i]=='-c' ){CONDITIONAL = as.logical(cmd[i+1])}
}

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

#######################################################
###### process all states to and combine eQTLs  #######
#######################################################

# for each gene, consider all SNPs in 95% CI for at least 1 component, celltype & condition
# for each SNP sum |Z| over significant component, celltype & condition, & effect direction
# select SNP with the highest score
# remove component, celltype & condition, & effect direction that include this SNP
# repeat until everything has been removed

QTL_assoc=list()
FDR_assoc=list()
for (CELLTYPE__STATE in cellstates){
  CELLTYPE=gsub('(.*)__(.*)','\\1',CELLTYPE__STATE)
  STATE=gsub('(.*)__(.*)','\\2',CELLTYPE__STATE)

  cat('\n\n\n',CELLTYPE,'-',STATE,'\n\n')
  QTL_assoc[[CELLTYPE__STATE]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/SusieR_95pct_CredibleSets_lbfover3.txt.gz',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),show_col_types = FALSE))
  if(CONDITIONAL){
    FDR_assoc[[CELLTYPE__STATE]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm_zvalcond.txt',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),show_col_types = FALSE))
  }else{
    FDR_assoc[[CELLTYPE__STATE]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm.txt',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),show_col_types = FALSE))
  }
  if(any(class(QTL_assoc[[CELLTYPE__STATE]])=='try-error')){
    cat('error assoc',CELLTYPE__STATE)
    QTL_assoc[[CELLTYPE__STATE]]=NULL
  }
  if(any(class(FDR_assoc[[CELLTYPE__STATE]])=='try-error')){
    cat('error FDR',CELLTYPE__STATE)
    FDR_assoc[[CELLTYPE__STATE]]=NULL
  }

}
QTL_assoc=rbindlist(QTL_assoc)
FDR_assoc=rbindlist(FDR_assoc,idcol='cellstate')


# for each SNP sum |Z| over significant component, celltype & condition, & effect direction
# select SNP with the highest score
# remove component, celltype & condition, & effect direction that include this SNP
# repeat until everything has been removed

QTL_assoc=merge(QTL_assoc, FDR_assoc[,.(cellstate,gene, top_component,FDR,dist_FDR)],by=c('cellstate','gene','top_component'))
QTL_assoc[,sign:=sign(zval)]
fwrite(QTL_assoc,sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT),sep='\t')
#QTL_assoc=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))


correct_factor=ifelse(CORRECT_CELLTYPE,length(cellstates),1)
bestSNP_list=NULL
QTL_assoc_remaining=QTL_assoc[FDR<FDR_TH/correct_factor & converged==TRUE,]
FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)
CONDITIONAL_char=ifelse(CONDITIONAL,'_zvalcond','')
CORRECT_CELLTYPE_char=ifelse(CORRECT_CELLTYPE,'_Nctadj','')

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
explained_components=merge(QTL_assoc[FDR<FDR_TH/correct_factor,],bestSNP_list,by=c('gene','snps'))
explained_components[,length(unique(snps)),by=gene][,.N,keyby=V1][,.(V1,N,N/sum(N))]
explained_components[,length(unique(cellstate)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
explained_components[,length(unique(celltype)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
explained_components[,length(unique(celltype)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]

fwrite(explained_components,file=sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))


################################################################################
####      merge results of eQTL mapping by cell types & joint eQTL mapping   ###
################################################################################


QTL_combined=list()
for (CHR in 1:22){
  cat(CHR,'\n')
  QTL_combined[[CHR]]=read_tsv(sprintf('%s/%s/dist_%s/FineMapping/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),show_col_types = FALSE)
}
QTL_combined=rbindlist(QTL_combined,idcol='CHR')
FDR_combined=try(read_tsv(file=sprintf('%s/%s/dist_%s/FDR_estimates_SusieR_perm_global.txt',OUT_DIR,RUN_NAME,CIS_DIST_TEXT),show_col_types = FALSE))
FDR_combined=as.data.table(FDR_combined)

QTL_assoc_combined=merge(QTL_combined[lbf_top_component>3 & minCSlevel_top_component<.95 & converged==TRUE,],FDR_combined,by=c('gene','top_component','converged'),suffix=c('.snp','.component'))
QTL_assoc_combined=QTL_assoc_combined[FDR<FDR_TH & converged==TRUE,]
QTL_assoc_combined[,cellstate:='combined']
QTL_assoc_combined[,zval:=combZ.snp]

# getPIP=function(x){N=table(x);rep(diff(c(0,sort(unique(x))))/N,N)}
#  QTL_assoc_combined[,pip_top_component:=getPIP(minCSlevel_top_component),by=.(gene,top_component)]

cols_to_keep=c('cellstate','gene','snps','top_component','pip_top_component','rank_top_component','minCSlevel_top_component','lbf_top_component','zval','FDR')

bestSNP_list=NULL
QTL_assoc_remaining0=rbind(QTL_assoc[FDR<FDR_TH/correct_factor & converged==TRUE,mget(cols_to_keep)],QTL_assoc_combined[FDR<FDR_TH & converged==TRUE,mget(cols_to_keep)])
QTL_assoc_remaining=QTL_assoc_remaining0
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
explained_components=merge(QTL_assoc_remaining0,bestSNP_list,by=c('gene','snps'))
explained_components[,celltype:=gsub('(.*)__(.*)','\\1',cellstate)]
explained_components[,state:=gsub('(.*)__(.*)','\\2',cellstate)]

explained_components[,length(unique(snps)),by=gene][,.N,keyby=V1][,.(V1,N,N/sum(N))]
explained_components[,length(unique(cellstate)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
explained_components[,length(unique(celltype)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
fwrite(explained_components,file=sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_and_joint_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
#fwrite(QTL_combined[minCSlevel_top_component<.95 & lbf_top_component>3,],file=sprintf('%s/%s/dist_%s/FineMapping/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))


#   QTL_combined=QTL_assoc[,.(bestP=min(pvalue),combZ=stouffer.method(pvalue,beta)),keyby=.(gene,snps)]
#   QTL_combined[,combP:=2*pnorm(abs(combZ),low=F)]
#   QTL_randSNP=QTL_combined[,.SD[sample(1:.N,1),],keyby=gene]
#   QTL_randSNP=merge(QTL_assoc,QTL_randSNP,by=c('gene','snps'))
#   fwrite(QTL_randSNP, sprintf('%s/%s/dist_%s/eQTL_randomSNP_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
#   QTL_combined=QTL_combined[order(gene,combP),]
#   fwrite(QTL_combined, sprintf('%s/%s/dist_%s/eQTL_combined_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
#   QTL_bestSNP=unique(QTL_combined,by="gene")
#   QTL_bestSNP=merge(QTL_assoc,QTL_bestSNP,by=c('gene','snps'))
#   fwrite(QTL_bestSNP, sprintf('%s/%s/dist_%s/eQTL_bestSNP_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
#
#   QTL_combined_perm=QTL_perm[,.(combZ=stouffer.method(pvalue,beta)),keyby=.(gene,snps)]
#   QTL_combined_perm[,combP:=2*pnorm(abs(combZ),low=F)]
#   fwrite(QTL_combined_perm, sprintf('%s/%s/dist_%s/eQTL_combined_chr%s_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
#   QTL_randSNP_perm=QTL_combined_perm[,.SD[sample(1:.N,1),],keyby=gene]
#   QTL_randSNP_perm=merge(QTL_perm,QTL_randSNP_perm,by=c('gene','snps'))
#   fwrite(QTL_randSNP_perm, sprintf('%s/%s/dist_%s/eQTL_randomSNP_chr%s_assoc_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
#   QTL_combined_perm=QTL_combined_perm[order(gene,combP),]
#   QTL_bestSNP_perm=unique(QTL_combined_perm,by="gene")
#   QTL_bestSNP_perm=merge(QTL_perm,QTL_bestSNP_perm,by=c('gene','snps'))
#   fwrite(QTL_bestSNP_perm, sprintf('%s/%s/dist_%s/eQTL_bestSNP_chr%s_assoc_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
# # write 4 files per chromosomes (best SNP from observed data, best SNP on permuted data, random SNP from observed/permuted data)
# rm(QTL_assoc,QTL_perm,QTL_bestSNP,QTL_randSNP,QTL_bestSNP_perm,QTL_randSNP_perm,QTL_combined_perm)
# gc()
#
# Susie.obj=readRDS(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/B__COV/FineMapping_100kb/SuSiE/SuSiE_output_ENSG00000213512.RDS",EVO_IMMUNO_POP_ZEUS))
# fwrite(QTL_combined[minCSlevel_top_component<.95 & lbf_top_component>3,],file=sprintf('%s/%s/dist_%s/FineMapping/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
