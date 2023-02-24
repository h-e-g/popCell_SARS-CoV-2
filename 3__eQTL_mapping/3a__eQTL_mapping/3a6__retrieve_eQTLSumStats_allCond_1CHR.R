################################################################################
################################################################################
# File name: 3a6__retrieve_eQTLSumStats_allCond_1CHR.R
# Author: Y.A., M.R., M.ON.
################################################################################
###############################################################################
# Step: retrieve eQTL_SumStats from all (lineage | celltype) x condition
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

RUN_NAME="/lineage_condition___CellPropLineage_SVs_RUN1"
RUN_NAME_ANNOT="/celltype_condition___CellPropLineage_SVs_RUN1"
CIS_DIST=1e5
CHR=22
FDR_TH=0.01
ANNOT_char=''
RUN_ID='RUN1'

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run from which to obatin eQTLs
  if (cmd[i]=='--run_name_annotate' | cmd[i]=='-a' ){RUN_NAME_ANNOT = cmd[i+1]} # ID of the run to use for annotation
	if (cmd[i]=='--run_id' | cmd[i]=='-i' ){RUN_ID = cmd[i+1]}
  if (cmd[i]=='--name' | cmd[i]=='-n' ){ANNOT_char = cmd[i+1]} # names to use characterize stats addes in the name of the ouput file (eg. celltype, or lineage)
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # FDR to consider eQTLs
}

if(is.null(RUN_NAME_ANNOT)){
  RUN_NAME_ANNOT=RUN_NAME
}

EQTL_DIR = "3__eQTL_mapping/SumStats"
CIS_DIST_TEXT=paste0(CIS_DIST/1000,'kb')
FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)

dir.create(sprintf('%s/%s/',EQTL_DIR,RUN_ID))

dir.create(sprintf('%s/%s/allStats_byChr',EQTL_DIR,RUN_ID))
cellstates=dir(sprintf('%s/%s',EQTL_DIR,RUN_NAME_ANNOT),pattern='.*__(NS|COV|IAV)')

#######################################################################################
###### process all chromosomes to extract statistics from all independant eQTLs  ######
#######################################################################################

  cat('\n\n\n',CHR,'\n\n')
  QTL_assoc=list()
  QTL_perm=list()
  for(cs in cellstates){
      cat('\n',cs,'\n')
        QTL_assoc[[paste(cs,CHR)]]=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc.txt.gz',EQTL_DIR,RUN_NAME_ANNOT,cs,CHR),show_col_types = FALSE))
        if(any(class(QTL_assoc[[paste(cs,CHR)]])=='try-error')){
          cat('error',cs,CHR)
          QTL_assoc[[paste(cs,CHR)]]=NULL
        }else{
          QTL_assoc[[paste(cs,CHR)]]$cellstate=cs
          QTL_assoc[[paste(cs,CHR)]]=as.data.table(QTL_assoc[[paste(cs,CHR)]])[CisDist<CIS_DIST,]
          }
        QTL_perm[[paste(cs,CHR)]]=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc_permuted.txt.gz',eQTL_DIR,RUN_NAME_ANNOT,cs,CHR),show_col_types = FALSE))
        if(any(class(QTL_perm[[paste(cs,CHR)]])=='try-error')){
          cat('error',cs,CHR)
          QTL_perm[[paste(cs,CHR)]]=NULL
        }else{
          QTL_perm[[paste(cs,CHR)]]$cellstate=cs
          QTL_perm[[paste(cs,CHR)]]=as.data.table(QTL_perm[[paste(cs,CHR)]])[CisDist<CIS_DIST,]
          }
      }
  QTL_assoc=rbindlist(QTL_assoc)
  QTL_perm=rbindlist(QTL_perm)


  independent_eQTLs=fread(file=sprintf('%s/%s/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR.txt.gz',EQTL_DIR,RUN_ID,FDR_char))
  QTL_allcond=merge(QTL_assoc,unique(independent_eQTLs[,.(gene,snps)]),by=c("gene","snps"))
	QTL_allcond[,se:=beta/(1e-7+statistic)]
	QTL_allcond[,beta_low:=sign(beta)*pmax(0,abs(beta)-2*se)]
  QTL_allcond[,beta_high:=sign(beta)*pmax(0,abs(beta)+2*se)]
  fwrite(QTL_allcond, file=sprintf('%s/%s/allStats_byChr/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR_allStats%s_chr%s.txt.gz',EQTL_DIR,RUN_ID,FDR_char,ANNOT_char,CHR),sep='\t')
	rm(QTL_assoc,QTL_perm)
	gc()
q('no')
# Max 20G per CHR
# max: 2 hours
#qos fast
