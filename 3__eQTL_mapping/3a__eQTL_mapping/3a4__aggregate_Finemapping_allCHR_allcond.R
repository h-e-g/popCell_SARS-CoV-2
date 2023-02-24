################################################################################
################################################################################
# File name: 3a4__aggregate_FineMapping_allCHR_allcond.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggegrate all Fine Mapped eQTLs from 1 condition & compute FDR
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


EQTL_DIR = "3_eQTL_mapping"
DATA_DIR = sprintf("%s/data",EQTL_DIR)

RUN_NAME="/lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5
FDR_TH=0.01

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
}

CIS_DIST_TEXT=paste0(CIS_DIST/1000,'kb')
FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)

COMBINED_EQTL_DIR=sprintf('%s/%s/dist_%s',EQTL_DIR,RUN_NAME,CIS_DIST_TEXT)

dir.create(COMBINED_EQTL_DIR)
cellstates=dir(sprintf('%s/%s',EQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV)')

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

  FINEMAP_DIR=sprintf("%s/%s/%s/FineMapping_%s/",EQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT)
  cat('\n\n\n',CELLTYPE,'-',STATE,'\n\n')
  QTL_assoc[[CELLTYPE__STATE]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/SusieR_95pct_CredibleSets_lbfover3.txt.gz',FINEMAP_DIR),show_col_types = FALSE))
  FDR_assoc[[CELLTYPE__STATE]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm.txt',FINEMAP_DIR),show_col_types = FALSE))
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
fwrite(QTL_assoc,sprintf('%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',COMBINED_EQTL_DIR),sep='\t')


bestSNP_list=NULL
QTL_assoc_remaining=QTL_assoc[FDR<FDR_TH & converged==TRUE,]
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
explained_components=merge(QTL_assoc[FDR<FDR_TH,],bestSNP_list,by=c('gene','snps'))
fwrite(explained_components,file=sprintf('%s/independent_eQTLs_allcond_%spctFDR.txt.gz',COMBINED_EQTL_DIR,FDR_char))

#### controls (check distribution of number of eQTL/per gene)
# explained_components[,length(unique(snps)),by=gene][,.N,keyby=V1][,.(V1,N,N/sum(N))]
#### controls (check distribution of number of cellstate per eQTL)
# explained_components[,length(unique(cellstate)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
#### controls (check distribution of number of celltype per eQTL)
# explained_components[,length(unique(celltype)),by=snps][,.N,keyby=V1][,.(V1,N,N/sum(N))]
