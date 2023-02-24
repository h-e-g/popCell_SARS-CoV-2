################################################################################
################################################################################
# File name: 3a7__aggregate_eQTLSumStats_allCond_allCHR.R
# Author: Y.A., M.R., M.ON.
################################################################################
###############################################################################
# Step: combine eQTL_SumStats from step 3a6 across all CHR
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
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # FDR to consider eQTLs
}

if(is.null(RUN_NAME_ANNOT)){
  RUN_NAME_ANNOT=RUN_NAME
}

FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)

tic('get gene_name for eQTL annotation')
feature_toUse=fread(sprintf('%s/data/genes_to_use.tsv',DAT_EQTL_DIR),header=F)$V1

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


tic('read fine mapping results')
QTL_fineMapped=list()
for (CHR in 1:22){
  cat('chr',CHR,'\n')
  QTL_fineMapped[[CHR]]=read_tsv(sprintf('%s/%s/allStats_byChr/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR_allStats%s_chr%s.txt.gz',DAT_EQTL_DIR,RUN_ID,FDR_char,ANNOT_char,CHR),show_col_types = FALSE)
  }
QTL_fineMapped=rbindlist(QTL_fineMapped,idcol='CHR')
QTL_fineMapped[,state:=gsub('(.*)__(.*)','\\2',cellstate)]

QTL_fineMapped=merge(QTL_fineMapped,Feature_annot[,.(gene=gene_id,gene_name)],by='gene')

fwrite(QTL_fineMapped,file=sprintf('%s/%s/independent_eQTLs_allcond_celltype_and_lineageLevel_%spctFDR_allStats%s_allCHR_FDR.txt.gz',DAT_EQTL_DIR,RUN_ID,FDR_char,ANNOT_char),sep='\t')
# toc()
