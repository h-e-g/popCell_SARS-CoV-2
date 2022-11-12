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
CELLTYPE__STATE='B__COV'
#CHR=22

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--cellstate' | cmd[i]=='-t' ){CELLTYPE__STATE = cmd[i+1]} # CELLTYPE__STATE to test
  # if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test

}

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
#cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1
myCELLTYPE=gsub('(.*)__(.*)','\\1',CELLTYPE__STATE)
mySTATE=gsub('(.*)__(.*)','\\2',CELLTYPE__STATE)

RUN_NAME=gsub('/','',RUN_NAME)
CELLTYPE=gsub('(lineage|celltype)_(condition|activation)((_logFC)?)__(.*)_([0-9]+)$','\\1',RUN_NAME)
STATE=gsub('(lineage|celltype)_(condition|activation)((_logFC)?)__(.*)_([0-9]+)$','\\2',RUN_NAME)
LOGFC=gsub('(lineage|celltype)_(condition|activation)((_logFC)?)__(.*)_([0-9]+)$','\\3',RUN_NAME)
COVSET=gsub('(lineage|celltype)_(condition|activation)((_logFC)?)__(.*)_([0-9]+)$','\\5',RUN_NAME)
COVSET=gsub('^(logFC_)?_','',COVSET)
DATE=gsub('(lineage|celltype)_(condition|activation)((_logFC)?)__(.*)_([0-9]+)$','\\6',RUN_NAME)

cat(CELLTYPE,STATE,LOGFC,COVSET,DATE,'\n',sprintf('%s_%s%s__%s',CELLTYPE,STATE,LOGFC,COVSET))

  cat('\n\n\n',myCELLTYPE,'-',mySTATE,'\n\n')

############################################################
###### process all chromosomes to identify best SNPs  ######
############################################################
QTL_assoc=list()
QTL_perm=list()

for( CHR in 1:22){
  cat('\n\n\n',CHR,'\n\n')
  #QTL_assoc[[paste(clust,CHR)]]=try(fread(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc.txt.gz',OUT_DIR,RUN_NAME,clust,CHR),sep='\t'))
  QTL_assoc[[CHR]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,CHR),show_col_types = FALSE))
  if(any(class(QTL_assoc[[CHR]])=='try-error')){
          cat('error',CELLTYPE__STATE,CHR)
          QTL_assoc[[CHR]]=NULL
        }else{
          QTL_assoc[[CHR]]$cellstate=CELLTYPE__STATE
          QTL_assoc[[CHR]]=as.data.table(QTL_assoc[[CHR]])[CisDist<CIS_DIST,]
          ## Quick correction of rank_top_component, minCSlevel_top_component and addition of pip_top_component
          # not needed when using the newest version of 9d1
          N=length(unique(QTL_assoc[[CHR]][,gene]))
          i=0
          cat('processing',N,'genes:')
          for (GENE in unique(QTL_assoc[[CHR]][,gene])[1:N]){
            i=i+1
            cat(i,'')
            susie.obj=try(readRDS(file=sprintf('%s/%s/%s/FineMapping_%s/SuSiE/SuSiE_output_%s.RDS',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,GENE)))
            if(class(susie.obj)!='try-error'){
              rank_mat=apply(-susie.obj$alpha,1,rank,ties.method='random')
              sorted_alpha_mat=apply(susie.obj$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
              alpha_mat=sorted_alpha_mat
              for (j in 1:ncol(sorted_alpha_mat)){
                alpha_mat[,j]=sorted_alpha_mat[rank_mat[,j],j]
              }
              QTL_assoc[[CHR]][gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_assoc[[CHR]][gene==GENE,.(1:.N,top_component)])]]
              QTL_assoc[[CHR]][gene==GENE,minCSlevel_top_component:=alpha_mat[as.matrix(QTL_assoc[[CHR]][gene==GENE,.(1:.N,top_component)])]]
              QTL_assoc[[CHR]][gene==GENE,pip_top_component:=susie.obj$alpha[as.matrix(QTL_assoc[[CHR]][gene==GENE,.(top_component,1:.N)])]]
              }else{cat('failed ')}
            }
            ## END Quick correction
          }
  QTL_perm[[CHR]]=try(read_tsv(file=sprintf('%s/%s/%s/FineMapping_%s/eQTL_FineMapped_chr%s_permuted.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,CHR),show_col_types = FALSE))
        if(any(class(QTL_perm[[CHR]])=='try-error')){
          cat('error',CELLTYPE__STATE,CHR)
          QTL_perm[[CHR]]=NULL
        }else{
          QTL_perm[[CHR]]$cellstate=CELLTYPE__STATE
          QTL_perm[[CHR]]=as.data.table(QTL_perm[[CHR]])[CisDist<CIS_DIST,]
          ## Quick correction of rank_top_component, minCSlevel_top_component and addition of pip_top_component
          # not needed when using the newest version of 9d1
          # N=length(unique(QTL_perm[[CHR]][,gene]))
          # i=0
          # for (GENE in unique(QTL_perm[[CHR]][,gene])[1:N]){
          #   i=i+1
          #   cat(i,'')
          #   susie.perm=try(readRDS(file=sprintf('%s/%s/%s/FineMapping_%s/SuSiE/SuSiE_output_%s.RDS',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,GENE)))
          #   if(class(susie.perm)!='try-error'){
          #     rank_mat=apply(-susie.perm$alpha,1,rank,ties.method='random')
          #     sorted_alpha_mat=apply(susie.perm$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
          #     alpha_mat=sorted_alpha_mat
          #     for (j in 1:ncol(sorted_alpha_mat)){
          #       alpha_mat[,j]=sorted_alpha_mat[rank_mat[,j],j]
          #     }
          #     QTL_perm[[CHR]][gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_perm[[CHR]][gene==GENE,.(1:.N,top_component)])]]
          #     QTL_perm[[CHR]][gene==GENE,minCSlevel_top_component:=alpha_mat[as.matrix(QTL_perm[[CHR]][gene==GENE,.(1:.N,top_component)])]]
          #     QTL_perm[[CHR]][gene==GENE,pip_top_component:=susie.perm$alpha[as.matrix(QTL_perm[[CHR]][gene==GENE,.(top_component,1:.N)])]]
          #     }else{cat('failed ')}
          #   }
            ## END Quick correction
          }
      }
  QTL_assoc=rbindlist(QTL_assoc)
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
  fwrite(QTL_zval[perm==0,], file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm.txt',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),sep='\t')

  minZ_th=QTL_zval[perm==0 & FDR<.05,min(abs(zval))]
  fwrite(QTL_assoc[minCSlevel_top_component<.95 & lbf_top_component>3 & pip_top_component>0.01,], file=sprintf('%s/%s/%s/FineMapping_%s/SusieR_95pct_CredibleSets_lbfover3.txt.gz',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),sep='\t')
  toc()

    tic('compare observed and permuted to estimate FDR - lbf')
    # move non converged to the end for observed and permuted to minimize their impact
    TP=unique(QTL_assoc[rank_top_component==1,.(gene,top_component,lbf_top_component,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(lbf_top_component))]
    TP[,perm:=0]
    FP=unique(QTL_perm[rank_top_component==1,.(gene,top_component,lbf_top_component,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(lbf_top_component))]
    FP[,perm:=1]
    QTL_lbf=rbind(TP,FP)[order(ifelse(converged,0,1),-abs(lbf_top_component)),]
    QTL_lbf$Nb_FP=cumsum(QTL_lbf$perm>0)/length(setdiff(QTL_lbf$perm,0))
    QTL_lbf$Nb_Pos=cumsum(QTL_lbf$perm==0)
    QTL_lbf$FDR=rev(cummin(rev(QTL_lbf$Nb_FP/QTL_lbf$Nb_Pos)))
    QTL_lbf[,dist_FDR:=CIS_DIST_TEXT]
    cat(QTL_lbf[perm==0 & FDR<.05,.N], 'eQTLs at 5% FDR', CIS_DIST/1000, 'kb')
    fwrite(QTL_lbf[perm==0,], file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm_lbf.txt',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),sep='\t')
  toc()


    tic('compare observed and permuted to estimate FDR - zvalconditional')
    # move non converged to the end for observed and permuted to minimize their impact
    TP=unique(QTL_assoc[rank_top_component==1,.(gene,top_component,t_conditional,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(t_conditional))]
    TP[,perm:=0]
    FP=unique(QTL_perm[rank_top_component==1,.(gene,top_component,t_conditional,converged)],by=c('gene','top_component'))[order(ifelse(converged,0,1),-abs(t_conditional))]
    FP[,perm:=1]
    QTL_zvalcond=rbind(TP,FP)[order(ifelse(converged,0,1),-abs(t_conditional)),]
    QTL_zvalcond$Nb_FP=cumsum(QTL_zvalcond$perm>0)/length(setdiff(QTL_zvalcond$perm,0))
    QTL_zvalcond$Nb_Pos=cumsum(QTL_zvalcond$perm==0)
    QTL_zvalcond$FDR=rev(cummin(rev(QTL_zvalcond$Nb_FP/QTL_zvalcond$Nb_Pos)))
    QTL_zvalcond[,dist_FDR:=CIS_DIST_TEXT]
    cat(QTL_zvalcond[perm==0 & FDR<.05,.N], 'eQTLs at 5% FDR', CIS_DIST/1000, 'kb')
    fwrite(QTL_zvalcond[perm==0,], file=sprintf('%s/%s/%s/FineMapping_%s/FDR_estimates_SusieR_perm_zvalcond.txt',OUT_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT),sep='\t')
    toc()

q('no')

# #### set up computation of LD adjusted Zvalues
# y=QTL_assoc[minCSlevel_top_component<.95 & lbf_top_component>3 & pip_top_component>0.01,]
# source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))
#
# Map=getMap(annotate=TRUE)
# SNP_list=getSNP(y[rank_top_component==1,unique(snps)],Map=Map)
#
# ##### define IID to use:
# MinCell_perCOND=500
# keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
# ##### load expression
# EXPRESSION_FILE=sprintf('%s/2_population_differences/BatchAdjusted_logCPM_125libs__per_%s_%s_annotated.tsv.gz',DATA_DIR,CELLTYPE,STATE)
# Expression=fread(file=EXPRESSION_FILE)
# Expression=Expression[,-c('ncells','Age','Gender')]
# # remove inds with < 500 cells
# Expression=Expression[IID%chin%keptIID,]
# for (GENE in unique(y$gene)){
#   # GENE="ENSG00000164308"
#   eQTL_snps=y[rank_top_component==1 & gene==GENE & cellstate==CELLTYPE__STATE,unique(snps)]
#   Expression_gene=merge(Expression[ID==GENE & celltype==myCELLTYPE & state==mySTATE,],dcast(SNP_list[ID%chin%eQTL_snps,],IID~ID,value.var='Number_of_ALT_alelle') ,by='IID')
#   set.seed(123)
#   rankTransform = function(x){
#       percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
#       mean_level=mean(x,na.rm=TRUE)
#       sd_level=sd(x,na.rm=TRUE)
#       qnorm(percentile,mean_level,sd_level)
#       }
#       Expression_gene[,logCPM:=rankTransform(logCPM)]
#       Expression_gene[,lm()]
#
#       COV_DIR=sprintf("%s/2_population_differences/Covariates",DATA_DIR)
#       COV_RUN_NAME=sprintf('%s_%s%s__%s',CELLTYPE,STATE,LOGFC,COVSET)
#       Covariates=fread(sprintf('%s/%s/Covariates__%s_%s.tsv.gz',COV_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE))
#       Expression_gene_withCOV=merge(Expression_gene[,-c('POP','celltype','state','Symbol','ID')],Covariates,by='IID')
#       # Expression_gene_withCOV[,summary(lm(logCPM~rs2927608+.,data=.SD[,-c('IID',"rs34028025","rs34358677","rs56285155","rs9918162","rs73152163","rs80041669","rs62377802",  "rs76241969",  "rs113210801"))]))]
#
#   }
#   CELLTYPE=gsub('(lineage|celltype)_(condition|activation)(.*)___(.*)_([0-9]+)$','\\1',RUN_NAME)
#   STATE=gsub('(lineage|celltype)_(condition|activation)(.*)___(.*)_([0-9]+)$','\\2',RUN_NAME)
#   LOGFC=gsub('(lineage|celltype)_(condition|activation)(.*)___(.*)_([0-9]+)$','\\3',RUN_NAME)
#
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

# Susie.obj=readRDS(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/B__COV/FineMapping_100kb/SuSiE/SuSiE_output_ENSG00000213512.RDS",EVO_IMMUNO_POP_ZEUS))
# fwrite(QTL_combined[minCSlevel_top_component<.95 & lbf_top_component>3,],file=sprintf('%s/%s/dist_%s/FineMapping/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
# for each gene, consider all SNPs in 95% CI for at least 1 component, celltype & condition
# for each SNP sum |Z| over significant component, celltype & condition, & effect direction
# select SNP with the highest score
# remove component, celltype & condition, & effect direction that include this SNP
# repeat until everything has been removed
