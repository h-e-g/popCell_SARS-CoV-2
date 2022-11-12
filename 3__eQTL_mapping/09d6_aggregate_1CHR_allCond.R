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

RUN_NAME="/lineage_condition___CellPropLineage_SVs_240322"
RUN_NAME_ANNOT=NULL
CIS_DIST=1e5
CHR=22
FDR_TH=0.01
CORRECT_CELLTYPE=FALSE
CONDITIONAL=FALSE
ANNOT_char=''
USE_JOINT=FALSE
USE_CELLTYPE=TRUE

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run from which to obatin eQTLs
  if (cmd[i]=='--run_name_annotate' | cmd[i]=='-a' ){RUN_NAME_ANNOT = cmd[i+1]} # ID of the run to use for annotation
  if (cmd[i]=='--name' | cmd[i]=='-n' ){ANNOT_char = cmd[i+1]}
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--fdr' | cmd[i]=='-f' ){FDR_TH = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--correct_celltype' | cmd[i]=='-t' ){CORRECT_CELLTYPE = as.logical(cmd[i+1])}
  if (cmd[i]=='--conditional' | cmd[i]=='-o' ){CONDITIONAL = as.logical(cmd[i+1])}
  if (cmd[i]=='--add_joint_mapping' | cmd[i]=='-j' ){USE_JOINT = as.logical(cmd[i+1])}
  if (cmd[i]=='--add_celltype_mapping' | cmd[i]=='-j' ){USE_CELLTYPE = as.logical(cmd[i+1])}
}
if(is.null(RUN_NAME_ANNOT)){
  RUN_NAME_ANNOT=RUN_NAME
}

FDR_char=substr(format(FDR_TH,scientific=FALSE),4,100)
CONDITIONAL_char=ifelse(CONDITIONAL,'_zvalcond','')
CORRECT_CELLTYPE_char=ifelse(CORRECT_CELLTYPE,'_Nctadj','')
JOINT_char=ifelse(USE_JOINT,'_allcond_and_joint','_allcond')
CELLTYPE_char=ifelse(USE_CELLTYPE,'_celltype_and_lineageLevel','')

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
dir.create(sprintf('%s/%s/dist_%s/allStats_byChr',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME_ANNOT),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

fisher.method=function(pval){
  S=sum(-2*log(pval))
  pchisq(S,2*length(pval),low=F)
}

stouffer.method=function(pval,beta){
  zval=sign(beta)*abs(qnorm(pval/2))
  Z=sum(zval)/sqrt(length(zval))
  #2*pnorm(abs(S),low=F)
  Z
}

#######################################################################################
###### process all chromosomes to extract statistics from all independant eQTLs  ######
#######################################################################################

  cat('\n\n\n',CHR,'\n\n')
  QTL_assoc=list()
  QTL_perm=list()
  for(cs in cellstates){
      cat('\n',cs,'\n')
        #QTL_assoc[[paste(clust,CHR)]]=try(fread(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,clust,CHR),sep='\t'))
        QTL_assoc[[paste(cs,CHR)]]=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc.txt.gz',eQTL_DIR,RUN_NAME_ANNOT,cs,CHR),show_col_types = FALSE))
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

  #### reactivate first time running
  # QTL_combined=QTL_assoc[,.(bestP=min(pvalue),combZ=stouffer.method(pvalue,beta)),keyby=.(gene,snps)]
  # QTL_combined[,combP:=2*pnorm(abs(combZ),low=F)]
  #
  # QTL_randSNP=QTL_combined[,.SD[sample(1:.N,1),],keyby=gene]
  # QTL_randSNP=merge(QTL_assoc,QTL_randSNP,by=c('gene','snps'))
  # fwrite(QTL_randSNP, sprintf('%s/%s/dist_%s/eQTL_randomSNP_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
  #
  # QTL_combined=QTL_combined[order(gene,combP),]
  # fwrite(QTL_combined, sprintf('%s/%s/dist_%s/eQTL_combined_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
  #
  # QTL_bestSNP=unique(QTL_combined,by="gene")
  # QTL_bestSNP=merge(QTL_assoc,QTL_bestSNP,by=c('gene','snps'))
  # fwrite(QTL_bestSNP, sprintf('%s/%s/dist_%s/eQTL_bestSNP_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
  #
  # QTL_combined_perm=QTL_perm[,.(combZ=stouffer.method(pvalue,beta)),keyby=.(gene,snps)]
  # QTL_combined_perm[,combP:=2*pnorm(abs(combZ),low=F)]
  # fwrite(QTL_combined_perm, sprintf('%s/%s/dist_%s/eQTL_combined_chr%s_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
  #
  # QTL_randSNP_perm=QTL_combined_perm[,.SD[sample(1:.N,1),],keyby=gene]
  # QTL_randSNP_perm=merge(QTL_perm,QTL_randSNP_perm,by=c('gene','snps'))
  # fwrite(QTL_randSNP_perm, sprintf('%s/%s/dist_%s/eQTL_randomSNP_chr%s_assoc_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')
  #
  # QTL_combined_perm=QTL_combined_perm[order(gene,combP),]
  # QTL_bestSNP_perm=unique(QTL_combined_perm,by="gene")
  # QTL_bestSNP_perm=merge(QTL_perm,QTL_bestSNP_perm,by=c('gene','snps'))
  # fwrite(QTL_bestSNP_perm, sprintf('%s/%s/dist_%s/eQTL_bestSNP_chr%s_assoc_permuted.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,CHR),sep='\t')

  independent_eQTLs=fread(file=sprintf('%s/%s/dist_%s/independent_eQTLs%s%s_%spctFDR%s%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char))
  QTL_allcond=merge(QTL_assoc,unique(independent_eQTLs[,.(gene,snps)]),by=c("gene","snps"))
  fwrite(QTL_allcond, file=sprintf('%s/%s/dist_%s/allStats_byChr/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char,CHR),sep='\t')
  QTL_allcond_perm=merge(QTL_perm,unique(independent_eQTLs[,.(gene,snps)]),by=c("gene","snps"))
  fwrite(QTL_allcond_perm, file=sprintf('%s/%s/dist_%s/allStats_byChr/independent_eQTLs%s%s_%spctFDR%s%s_allStats%s_permuted_chr%s.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT,JOINT_char,CELLTYPE_char,FDR_char,CONDITIONAL_char,CORRECT_CELLTYPE_char,ANNOT_char,CHR),sep='\t')

  #QTL_allcond=merge(QTL_combined,unique(independent_eQTLs[,.(gene,snps)]),by=c("gene","snps"))

# write 4 files per chromosomes (best SNP from observed data, best SNP on permuted data, random SNP from observed/permuted data)
rm(QTL_assoc,QTL_perm,QTL_bestSNP,QTL_randSNP,QTL_bestSNP_perm,QTL_randSNP_perm,QTL_combined_perm)
gc()
q('no')
# Max 20G per CHR
# max: 2 hours
#qos fast
