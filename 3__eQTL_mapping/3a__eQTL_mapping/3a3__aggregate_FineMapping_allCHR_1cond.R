################################################################################
################################################################################
# File name: 3a3__aggregate_FineMapping_allCHR_1cond.R
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

# set defaults values
RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5
CHR=22
CELLTYPE='lineage'
COV_RUN_NAME="lineage_condition___CellPropLineage_SVs"

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
  if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use to define pseudoBulk. should match what was used for matrix eQTL
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COVAR_DIR/Covariates). should match what was used for matrix eQTL.
  if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run (folder where to find matrix eQTL output)
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--cellstate' | cmd[i]=='-t' ){CELLTYPE__STATE = cmd[i+1]} # CELLTYPE__STATE to test
}

myCELLTYPE=gsub('(.*)__(.*)','\\1',CELLTYPE__STATE)
mySTATE=gsub('(.*)__(.*)','\\2',CELLTYPE__STATE)
CIS_DIST_TEXT=paste0(CIS_DIST/1000,'kb')

EQTL_DIR=sprintf('3_eQTL_mapping/sumStats/%s/%s/',RUN_NAME,CELLTYPE__STATE)
FINEMAP_DIR=sprintf('%s/FineMapping_%s',EQTL_DIR,CIS_DIST_TEXT)

cat('\n\n\n',myCELLTYPE,'-',mySTATE,'\n\n')

############################################################
###### process all chromosomes to identify best SNPs  ######
############################################################

###### observed data
QTL_assoc=list()
for( CHR in 1:22){
  cat('\n\n\n',CHR,'\n\n')
  QTL_assoc[[CHR]]=try(read_tsv(file=sprintf('%s/eQTL_FineMapped_chr%s.txt.gz',FINEMAP_DIR,CHR),show_col_types = FALSE))
  if(any(class(QTL_assoc[[CHR]])=='try-error')){
    cat('error',CELLTYPE__STATE,CHR)
    QTL_assoc[[CHR]]=NULL
  }else{
    QTL_assoc[[CHR]]$cellstate=CELLTYPE__STATE
    QTL_assoc[[CHR]]=as.data.table(QTL_assoc[[CHR]])[CisDist<CIS_DIST,]
  }
}
QTL_assoc=rbindlist(QTL_assoc)

###### permuted data
QTL_perm=list()
for( CHR in 1:22){
  QTL_perm[[CHR]]=try(read_tsv(file=sprintf('%s/eQTL_FineMapped_chr%s_permuted.txt.gz',FINEMAP_DIR,CHR),show_col_types = FALSE))
  if(any(class(QTL_perm[[CHR]])=='try-error')){
    cat('error',CELLTYPE__STATE,CHR)
    QTL_perm[[CHR]]=NULL
  }else{
    QTL_perm[[CHR]]$cellstate=CELLTYPE__STATE
    QTL_perm[[CHR]]=as.data.table(QTL_perm[[CHR]])[CisDist<CIS_DIST,]
  }
}
QTL_perm=rbindlist(QTL_perm)

tic('compare observed and permuted to estimate FDR')
# move non converged to the end for observed and permuted to minimize their impact
TP=unique(QTL_assoc[rank_top_component==1,.(gene,top_component,zval,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(zval))]
TP[,perm:=0]
FP=unique(QTL_perm[rank_top_component==1,.(gene,top_component,zval,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(zval))]
FP[,perm:=1]
QTL_zval=rbind(TP,FP)[order(ifelse(converged,0,1),-abs(zval)),]
QTL_zval$Nb_FP=cumsum(QTL_zval$perm>0)/length(setdiff(QTL_zval$perm,0))
QTL_zval$Nb_Pos=cumsum(QTL_zval$perm==0)
QTL_zval$FDR=rev(cummin(rev(QTL_zval$Nb_FP/QTL_zval$Nb_Pos)))
QTL_zval[,dist_FDR:=CIS_DIST_TEXT]
cat(QTL_zval[perm==0 & FDR<.05,.N], 'eQTLs at 5% FDR', CIS_DIST/1000, 'kb')
fwrite(QTL_zval[perm==0,], file=sprintf('%s/FDR_estimates_SusieR_perm.txt',FINEMAP_DIR),sep='\t')

minZ_th=QTL_zval[perm==0 & FDR<.05,min(abs(zval))]
fwrite(QTL_assoc[minCSlevel_top_component<.95 & lbf_top_component>3 & pip_top_component>0.01,], file=sprintf('%s/SusieR_95pct_CredibleSets_lbfover3.txt.gz',FINEMAP_DIR),sep='\t')
toc()

q('no')
