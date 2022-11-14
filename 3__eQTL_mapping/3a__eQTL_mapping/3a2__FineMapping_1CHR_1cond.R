options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer))
suppressMessages(require(susieR))

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
source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))

RUN_NAME="/lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5
CHR=22

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--cellstate' | cmd[i]=='-t' ){CELLTYPE__STATE = cmd[i+1]} # CELLTYPE__STATE to test
}

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
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

GET_LOGFC=ifelse(LOGFC=='_logFC',TRUE,FALSE)

cat('\n\n\n',CHR,'-',myCELLTYPE,'-',mySTATE,'\n\n')
QTL_assoc=list()
QTL_perm=list()
QTL_assoc=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CHR),show_col_types = FALSE))
if(any(class(QTL_assoc)=='try-error')){
  cat('error',CELLTYPE__STATE,CHR)
  QTL_assoc=NULL
}else{
  QTL_assoc$cellstate=CELLTYPE__STATE
  QTL_assoc=as.data.table(QTL_assoc)[CisDist<CIS_DIST,]
  }

QTL_perm=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc_permuted.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CHR),show_col_types = FALSE))
if(any(class(QTL_perm)=='try-error')){
  cat('error',CELLTYPE__STATE,CHR)
  QTL_perm=NULL
  }else{
    QTL_perm$cellstate=CELLTYPE__STATE
    QTL_perm=as.data.table(QTL_perm)[CisDist<CIS_DIST,]
    }

QTL_assoc[,zval:=sign(beta)*abs(qnorm(pvalue/2))]
QTL_perm[,zval:=sign(beta)*abs(qnorm(pvalue/2))]

susie.method=function(zval,R){
#  zval=sign(beta)*abs(qnorm(pval/2))
  fitted_rss <- susie_rss(zval, R, L = 10, niter=200)
  fitted_rss
}

##### define IID to use:
MinCell_perCOND=500
filtIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/filteredIID_lessThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]

##### load genotypes
Genotypes=queryRange(CHR,0,300e6)
Genotypes=melt(Genotypes,id.vars=c("ID"),variable.name='IID')[IID%in%keptIID,]
setkey(Genotypes,ID,IID)
Genotypes[,value:=as.numeric(value)]

dir.create(sprintf('%s/%s/%s/FineMapping_%s/',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT))
dir.create(sprintf('%s/%s/%s/FineMapping_%s/SuSiE/',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT))

##### load expression
EXPRESSION_FILE=sprintf('%s/2_population_differences/BatchAdjusted_logCPM_125libs__per_%s_%s_annotated.tsv.gz',DATA_DIR,CELLTYPE,STATE)
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]
# remove inds with < 500 cells
Expression=Expression[IID%chin%keptIID,]

tic('select target cell type, state, individuals and features')
# select target cell type and state
if(GET_LOGFC==TRUE){
  Expression_NS=Expression[celltype==gsub('.INFECTED','',myCELLTYPE) & state=="NS",]
	Expression_NS[,celltype:=myCELLTYPE]
}
Expression=Expression[celltype==myCELLTYPE & state==mySTATE,]
# obtain logFC
if(GET_LOGFC==TRUE){
  Expression=merge(Expression,Expression_NS,by=c('IID','celltype','ID','Symbol','POP'),suffix=c('','.NS'))
	Expression[,logCPM:=logCPM-logCPM.NS]
}
toc()

set.seed(123)
rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
    }

Expression[,logCPM:=rankTransform(logCPM),by=.(ID,Symbol,celltype,state)]

####### load covariates
COV_DIR=sprintf("%s/2_population_differences/Covariates",DATA_DIR)
COV_RUN_NAME=sprintf('%s_%s%s__%s',CELLTYPE,STATE,LOGFC,COVSET)
Covariates=fread(sprintf('%s/%s/Covariates__%s_%s.tsv.gz',COV_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE))



################################### run susie on all genes
N=length(unique(QTL_assoc[,gene]))
i=0

for (GENE in unique(QTL_assoc[,gene])[1:N]){
  i=i+1
  #GENE=unique(QTL_assoc[,gene])[i]
  tic(paste(i,'/',N,':',GENE,'\n'))

  variants=QTL_assoc[gene==GENE,snps]
  if(length(variants)>50){
    Geno_mat=dcast(Genotypes[variants,],IID~ID)
    Geno_res=apply(as.matrix(Geno_mat[,-'IID']),2,function(x){lm(x~substr(Geno_mat$IID,1,3))$res})
    R=cor(Geno_res)
      # fine map eQTL for observed data
    susie.obj=susie_rss(QTL_assoc[gene==GENE,zval],R[variants,variants],coverage=.8)
    saveRDS(susie.obj,file=sprintf('%s/%s/%s/FineMapping_%s/SuSiE/SuSiE_output_%s.RDS',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,GENE))
    QTL_assoc[gene==GENE,top_component:=apply(susie.obj$alpha,2,which.max)]
    rank_mat=apply(-susie.obj$alpha,1,rank,ties.method='random')
    sorted_alpha_mat=apply(susie.obj$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
    alpha_mat=sorted_alpha_mat
    for (j in 1:ncol(sorted_alpha_mat)){
      alpha_mat[,j]=sorted_alpha_mat[rank_mat[,j],j]
    }
    #alpha_mat=apply(susie.obj$alpha,1,function(x){c(-cumsum(sort(-x))[rank(-x)]})
    QTL_assoc[gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_assoc[gene==GENE,.(1:.N,top_component)])]]
    QTL_assoc[gene==GENE,minCSlevel_top_component:=alpha_mat[as.matrix(QTL_assoc[gene==GENE,.(1:.N,top_component)])]]
    QTL_assoc[gene==GENE,pip_top_component:=susie.obj$alpha[as.matrix(QTL_assoc[gene==GENE,.(top_component,1:.N)])]]
    QTL_assoc[gene==GENE,lbf_top_component:=susie.obj$lbf[top_component]]
    #QTL_assoc[gene==GENE,mu_top_component:=susie.obj$mu[top_component]]
    #QTL_assoc[gene==GENE,mu_se_top_component:=sqrt(susie.obj$mu2[top_component]-susie.obj$mu[top_component]^2)]
    QTL_assoc[gene==GENE,converged:=rep(susie.obj$converged,.N)]

    ################################### run linear model to get conditional Z-values for each component

    eQTL_top=QTL_assoc[gene==GENE & rank_top_component==1,]
    eQTL_snps=eQTL_top[order(-abs(zval)),.(snps,top_component)]
    SNP_list=dcast(Genotypes[ID%chin%eQTL_snps$snps,],IID~ID,value.var='value')
    Expression_gene=merge(Expression[ID==GENE,],SNP_list,by='IID')
    Expression_gene_withCOV=merge(Expression_gene[,-c('POP','celltype','state','Symbol','ID')],Covariates,by='IID')
     adj.set=c()
     coeff=list()
     cov.set=setdiff(colnames(Covariates),'IID')
     for (j in 1:length(eQTL_snps$snps)){
         adj.set=c(adj.set,eQTL_snps$snps[j-1])
         target.snp=eQTL_snps$snps[j]
         if(sd(Expression_gene_withCOV[,get(target.snp)])!=0){
           snp_coeff=Expression_gene_withCOV[,summary(lm(logCPM~.,data=.SD[,mget(c('logCPM',target.snp,adj.set,cov.set))]))$coeff[target.snp,]]
           coeff[[j]]=data.table(order=j,snp=target.snp,stat=names(snp_coeff),value=snp_coeff)
         }else{
           coeff[[j]]=data.table(order=j,snp=target.snp,stat=c("Estimate","Std. Error","t value","Pr(>|t|)"),value=c(0,NA,0,1))
         }
     }
     coeff=dcast(rbindlist(coeff),order+snp~stat)
     coeff=coeff[match(eQTL_top$snps,coeff$snp),]
     setnames(coeff,c("Estimate","Std. Error","t value","Pr(>|t|)"),c('beta','se','tvalue','pvalue'),skip_absent=TRUE)
     mm=match(eQTL_top$snps,coeff$snp)
     QTL_assoc[gene==GENE & rank_top_component==1 ,beta_conditional:=coeff[mm,beta]]
     QTL_assoc[gene==GENE & rank_top_component==1 ,se_conditional:=coeff[mm,se]]
     QTL_assoc[gene==GENE & rank_top_component==1 ,t_conditional:=coeff[mm,tvalue]]

      # fine map eQTL for permuted data
    susie.perm=susie_rss(QTL_perm[gene==GENE,zval],R[variants,variants],coverage=.8)
    saveRDS(susie.perm,file=sprintf('%s/%s/%s/FineMapping_%s/SuSiE/SuSiE_output_%s_perm.RDS',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,GENE))
    QTL_perm[gene==GENE,top_component:=apply(susie.perm$alpha,2,which.max)]
    rank_mat=apply(-susie.perm$alpha,1,rank,ties.method='random')
    sorted_alpha_mat=apply(susie.obj$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
    alpha_mat=sorted_alpha_mat
    for (j in 1:ncol(sorted_alpha_mat)){
      alpha_mat[,j]=sorted_alpha_mat[rank_mat[,j],j]
    }
    # alpha_mat=apply(susie.perm$alpha,1,function(x){-cumsum(sort(-x))[rank(-x)]})
    QTL_perm[gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_perm[gene==GENE,.(1:.N,top_component)])]]
    QTL_perm[gene==GENE,minCSlevel_top_component:=alpha_mat[as.matrix(QTL_perm[gene==GENE,.(1:.N,top_component)])]]
    QTL_perm[gene==GENE,pip_top_component:=susie.perm$alpha[as.matrix(QTL_perm[gene==GENE,.(top_component,1:.N)])]]
    QTL_perm[gene==GENE,lbf_top_component:=susie.perm$lbf[top_component]]
    QTL_perm[gene==GENE,converged:=rep(susie.perm$converged,.N)]

    QTL_perm[gene==GENE,beta_conditional:=beta]
    QTL_perm[gene==GENE,se_conditional:=se]
    QTL_perm[gene==GENE,t_conditional:=statistic]

  }else{
    QTL_assoc[gene==GENE,rank_top_component:=NA]
    QTL_assoc[gene==GENE,minCSlevel_top_component:=NA]
    QTL_assoc[gene==GENE,pip_top_component:=NA]
    QTL_assoc[gene==GENE,lbf_top_component:=NA]
    QTL_assoc[gene==GENE,converged:=FALSE]

    QTL_assoc[gene==GENE,beta_conditional:=NA]
    QTL_assoc[gene==GENE,se_conditional:=NA]
    QTL_assoc[gene==GENE,t_conditional:=NA]

    QTL_perm[gene==GENE,rank_top_component:=NA]
    QTL_perm[gene==GENE,minCSlevel_top_component:=NA]
    QTL_perm[gene==GENE,pip_top_component:=NA]
    QTL_perm[gene==GENE,lbf_top_component:=NA]
    QTL_perm[gene==GENE,converged:=FALSE]

    QTL_perm[gene==GENE,beta_conditional:=NA]
    QTL_perm[gene==GENE,se_conditional:=NA]
    QTL_perm[gene==GENE,t_conditional:=NA]

  }
  toc()
}
converged=QTL_assoc[,.(converged=any(converged)),by=gene]
print(converged[,paste('FineMapping converged for',sum(converged),'/',length(converged),'genes')])

fwrite(QTL_assoc,file=sprintf('%s/%s/%s/FineMapping_%s/eQTL_FineMapped_chr%s.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,CHR),sep='\t')
fwrite(QTL_perm,file=sprintf('%s/%s/%s/FineMapping_%s/eQTL_FineMapped_chr%s_permuted.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CIS_DIST_TEXT,CHR),sep='\t')
