options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
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

RUN_NAME="/lineage_condition___CellPropLineage_SVs_240322"
RUN_NAME_ANNOT=NULL

CIS_DIST=1e5
FDR_TH=0.01
CORRECT_CELLTYPE=FALSE
CONDITIONAL=FALSE
ANNOT_char=''
USE_JOINT=FALSE
USE_CELLTYPE=FALSE

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--run_name_annotate' | cmd[i]=='-a' ){RUN_NAME_ANNOT = cmd[i+1]} # ID of the run to use for annotation
  if (cmd[i]=='--name' | cmd[i]=='-n' ){ANNOT_char = cmd[i+1]}
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--correct_celltype' | cmd[i]=='-t' ){CORRECT_CELLTYPE = as.logical(cmd[i+1])}
  if (cmd[i]=='--conditional' | cmd[i]=='-o' ){CONDITIONAL = as.logical(cmd[i+1])}
  if (cmd[i]=='--add_joint_mapping' | cmd[i]=='-j' ){USE_JOINT = as.logical(cmd[i+1])}
  if (cmd[i]=='--add_celltype_mapping' | cmd[i]=='-j' ){USE_CELLTYPE = as.logical(cmd[i+1])}
}

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
if(is.null(RUN_NAME_ANNOT)){
  RUN_NAME_ANNOT=RUN_NAME
}
FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)
CONDITIONAL_char=ifelse(CONDITIONAL,'_zvalcond','')
CORRECT_CELLTYPE_char=ifelse(CORRECT_CELLTYPE,'_Nctadj','')
JOINT_char=ifelse(USE_JOINT,'_allcond_and_joint','_allcond')
CELLTYPE_char=ifelse(USE_CELLTYPE,'_celltype_and_lineageLevel','')
dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

tic('read fine mapping results')
QTL_fineMapped=list()
QTL_fineMapped_perm=list()
for (CHR in 1:22){
  cat('chr',CHR,'\n')
  QTL_fineMapped[[CHR]]=read_tsv(sprintf('%s/%s/dist_%s/allStats_byChr/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char,CHR),show_col_types = FALSE)
  QTL_fineMapped_perm[[CHR]]=read_tsv(sprintf('%s/%s/dist_%s/allStats_byChr/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_permuted_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char,CHR),show_col_types = FALSE)
  }
QTL_fineMapped=rbindlist(QTL_fineMapped,idcol='CHR')
QTL_fineMapped_perm=rbindlist(QTL_fineMapped_perm,idcol='CHR')
QTL_fineMapped[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
QTL_fineMapped_perm[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
fwrite(QTL_fineMapped,file=sprintf('%s/%s/dist_%s/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_allCHR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char),sep='\t')
fwrite(QTL_fineMapped_perm,file=sprintf('%s/%s/dist_%s/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_permuted_allCHR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char),sep='\t')

toc()

tic('get gene_name and annotate')
gtf <- rtracklayer::import(sprintf("%s/single_cell/resources/references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
Feature_annot=as.data.table(gtf)[type=='gene',.(gene_id,gene_name,seqnames, start, end, strand,gene_type)]
Feature_annot[is.na(gene_name),gene_name:=gene_id]
Feature_annot[,seqnames:=as.character(seqnames)]
# for IAV_M and IAV_NS  we have two lines with the same genename (2 transcripts)
# we only keep the longest transcript
Feature_annot=Feature_annot[!(gene_id=="IAV_M" & end==784) & !(gene_id=="IAV_NS" & end==719),]
# we match the remaining features to feature_toUse
Feature_annot=Feature_annot[match(feature_toUse,gene_id),]

QTL_fineMapped=merge(QTL_fineMapped,Feature_annot[,.(gene=gene_id,gene_name)],by='gene')
QTL_fineMapped_perm=merge(QTL_fineMapped_perm,Feature_annot[,.(gene=gene_id,gene_name)],by='gene')
fwrite(QTL_fineMapped,file=sprintf('%s/%s/dist_%s/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_allCHR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char),sep='\t')
fwrite(QTL_fineMapped_perm,file=sprintf('%s/%s/dist_%s/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_permuted_allCHR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char),sep='\t')
toc()

tic('create beta and se matrices')
setnames(QTL_fineMapped,c('gene_id'),c('gene'),skip_absent=TRUE)
# observed best SNPs
QTL_fineMapped[,se:=beta/(1e-7+statistic)]
beta_mat=dcast(QTL_fineMapped,gene+snps~cellstate,value.var='beta',fill=0) # if beta is NA, force it to 0
se_mat=dcast(QTL_fineMapped,gene+snps~cellstate,value.var='se',fill=1e-7) # if se is NA, force it to epsilon>0
toc()

# add minimal amount of noise to avoid numerical issues
Smin=1e-7;
add.noise=function(Matrix,noise.sd,rownames){
	set.seed(1234)
	nr=nrow(Matrix)
	nc=ncol(Matrix)
	Matrix = as.matrix(Matrix) + matrix(rnorm(nr*nc,0,noise.sd),nr,nc)
	rownames(Matrix)=rownames
	Matrix
}

id.vars=c('gene','snps')

########################################################################
##################### 				RUNNING MASHR 			######################
########################################################################

library(mashr)
tic('prepare mash input and estimate null correlations')

data.noised = mash_set_data(add.noise(beta_mat[,!..id.vars], Smin, paste(beta_mat[,gene],beta_mat[,snps],sep='_')),
                              as.matrix(se_mat[,!..id.vars])+Smin)

Vhat = readRDS(file=sprintf('%s/%s/dist_%s/mash_null_correlation_effectsSizes.RDS',OUT_DIR,RUN_NAME_ANNOT,CIS_DIST_TEXT))
data.Vhat = mash_update_data(data.noised, V=Vhat)
toc()
m=readRDS(file=sprintf('%s/%s/dist_%s/mash_fitted_null_model.RDS',OUT_DIR,RUN_NAME_ANNOT,CIS_DIST_TEXT))

tic('running mash on the set of top SNPs, (observed data)')
m2 = mash(data.Vhat, g=get_fitted_g(m), fixg=TRUE)
lfsr_mashr=get_lfsr(m2)
lfsr_mashr=data.table(reshape2::melt(lfsr_mashr))
colnames(lfsr_mashr)=c('gene_snp','cellstate','lfsr')
lfsr_mashr[,gene_id:=gsub('(.*)_(.*)','\\1',gene_snp)]
lfsr_mashr[,snps:=gsub('(.*)_(.*)','\\2',gene_snp)]

beta_mashr=get_pm(m2)
beta_mashr=data.table(reshape2::melt(beta_mashr))
colnames(beta_mashr)=c('gene_snp','cellstate','beta_mash')
beta_mashr[,gene_id:=gsub('(.*)_(.*)','\\1',gene_snp)]
beta_mashr[,snps:=gsub('(.*)_(.*)','\\2',gene_snp)]

se_mashr=get_psd(m2)
se_mashr=data.table(reshape2::melt(se_mashr))
colnames(se_mashr)=c('gene_snp','cellstate','se_mash')
se_mashr[,gene_id:=gsub('(.*)_(.*)','\\1',gene_snp)]
se_mashr[,snps:=gsub('(.*)_(.*)','\\2',gene_snp)]

saveRDS(m2,file=sprintf('%s/%s/dist_%s/mash_estimated_params_observed_independent_eQTLs%s%s_%spctFDR%s%s_allStats%s.RDS',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char))
toc()

setnames(QTL_fineMapped,c('gene'),c('gene_id'),skip_absent=TRUE)
QTL_fineMapped=merge(QTL_fineMapped,lfsr_mashr[,.(gene_id,snps,cellstate,lfsr)],by=c('cellstate','gene_id','snps'),all.x=TRUE)
QTL_fineMapped=merge(QTL_fineMapped,beta_mashr[,.(gene_id,snps,cellstate,beta_mash)],by=c('cellstate','gene_id','snps'),all.x=TRUE)
QTL_fineMapped=merge(QTL_fineMapped,se_mashr[,.(gene_id,snps,cellstate,se_mash)],by=c('cellstate','gene_id','snps'),all.x=TRUE)
fwrite(QTL_fineMapped,file=sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_mashr_and_FDR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char),sep='\t')
