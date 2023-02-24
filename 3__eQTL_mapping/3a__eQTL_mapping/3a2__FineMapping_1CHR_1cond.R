################################################################################
################################################################################
# File name: 3a2__FineMapping_1CHR_1cond.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use SusieR to fineMap eQTLs in a 100kb window
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
source(sprintf("MISC/misc_plots.R"))
source(sprintf("MISC/querySNPs.R"))

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
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test.
  if (cmd[i]=='--cellstate' | cmd[i]=='-t' ){CELLTYPE__STATE = cmd[i+1]} # CELLTYPE__STATE to test
}

##### define IID & genes to use:
IID_to_use=fread('3__eQTL_mapping/data/IID_individual_to_use.tsv',sep='\t',header=F)[,V1]
feature_toUse=fread('3__eQTL_mapping/data/genes_to_use.tsv',header=F)$V1

myCELLTYPE=gsub('(.*)__(.*)','\\1',CELLTYPE__STATE)
mySTATE=gsub('(.*)__(.*)','\\2',CELLTYPE__STATE)
CIS_DIST_TEXT=paste0(CIS_DIST/1000,'kb')

EQTL_DIR=sprintf('3_eQTL_mapping/sumStats/%s/%s/',RUN_NAME,CELLTYPE__STATE)
FINEMAP_DIR=sprintf('%s/FineMapping_%s',EQTL_DIR,CIS_DIST_TEXT)

dir.create(FINEMAP_DIR)
dir.create(sprintf('%s/SuSiE/',FINEMAP_DIR))

cat('\n\n\n',CHR,'-',myCELLTYPE,'-',mySTATE,'\n\n')
#### read association statistics (observed)
QTL_assoc=list()
QTL_assoc=try(read_tsv(file=sprintf('%s/eQTL_ALL_chr%s_assoc.txt.gz',EQTL_DIR,CHR),show_col_types = FALSE))
if(any(class(QTL_assoc)=='try-error')){
  cat('error',CELLTYPE__STATE,CHR)
  QTL_assoc=NULL
}else{
  QTL_assoc$cellstate=CELLTYPE__STATE
  QTL_assoc=as.data.table(QTL_assoc)[CisDist<CIS_DIST,]
  QTL_assoc[,zval:=sign(beta)*abs(qnorm(pvalue/2))]
  }

#### read association statistics (permuted)
QTL_perm=list()
QTL_perm=try(read_tsv(file=sprintf('%s/%s/%s/eQTL_ALL_chr%s_assoc_permuted.txt.gz',eQTL_DIR,RUN_NAME,CELLTYPE__STATE,CHR),show_col_types = FALSE))
if(any(class(QTL_perm)=='try-error')){
  cat('error',CELLTYPE__STATE,CHR)
  QTL_perm=NULL
  }else{
    QTL_perm$cellstate=CELLTYPE__STATE
    QTL_perm=as.data.table(QTL_perm)[CisDist<CIS_DIST,]
    QTL_perm[,zval:=sign(beta)*abs(qnorm(pvalue/2))]
}


susie.method=function(zval,R){
  fitted_rss <- susie_rss(zval, R, L = 10, niter=200)
  fitted_rss
}


##### load genotypes
Genotypes=queryRange(CHR,0,300e6)
Genotypes=melt(Genotypes,id.vars=c("ID"),variable.name='IID')[IID%in%IID_to_use,]
setkey(Genotypes,ID,IID)
Genotypes[,value:=as.numeric(value)]


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
    saveRDS(susie.obj,file=sprintf('%s/SuSiE/SuSiE_output_%s.RDS',FINEMAP_DIR,GENE))
    # keep top component for each SNP
    QTL_assoc[gene==GENE,top_component:=apply(susie.obj$alpha,2,which.max)]
    ### compute rank of each SNP on top component and cumulative probability for each SNP (to define 95% CI for each component)
    # rank SNPs by decreasing alpha (PIP) in each component
    rank_mat=apply(-susie.obj$alpha,1,rank,ties.method='random')
    # for SNP with the k_th highest PIP, compute cumulative_PIP_k = sum_i=1..k { alpha_(i) }
    sorted_cumulative_PIP=apply(susie.obj$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
    # reorder to have cumulative_PIP associated to each SNP with the original order
    cumulative_PIP=sorted_cumulative_PIP
    for (j in 1:ncol(sorted_cumulative_PIP)){
      cumulative_PIP[,j]=sorted_cumulative_PIP[rank_mat[,j],j]
    }
    QTL_assoc[gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_assoc[gene==GENE,.(1:.N,top_component)])]]
    QTL_assoc[gene==GENE,minCSlevel_top_component:=cumulative_PIP_k[as.matrix(QTL_assoc[gene==GENE,.(1:.N,top_component)])]]
    QTL_assoc[gene==GENE,pip_top_component:=susie.obj$alpha[as.matrix(QTL_assoc[gene==GENE,.(top_component,1:.N)])]]
    QTL_assoc[gene==GENE,lbf_top_component:=susie.obj$lbf[top_component]]
    QTL_assoc[gene==GENE,converged:=rep(susie.obj$converged,.N)]

    # fine map eQTL for permuted data
    susie.perm=susie_rss(QTL_perm[gene==GENE,zval],R[variants,variants],coverage=.8)
    saveRDS(susie.perm,file=sprintf('%s/SuSiE/SuSiE_output_%s_perm.RDS',FINEMAP_DIR,GENE))
    QTL_perm[gene==GENE,top_component:=apply(susie.perm$alpha,2,which.max)]
    rank_mat=apply(-susie.perm$alpha,1,rank,ties.method='random')
    sorted_cumulative_PIP=apply(susie.obj$alpha,1,function(x){c(0,-cumsum(sort(-x))[-length(x)])})
    cumulative_PIP=sorted_cumulative_PIP
    for (j in 1:ncol(sorted_cumulative_PIP)){
      cumulative_PIP[,j]=sorted_cumulative_PIP[rank_mat[,j],j]
    }
    QTL_perm[gene==GENE,rank_top_component:=rank_mat[as.matrix(QTL_perm[gene==GENE,.(1:.N,top_component)])]]
    QTL_perm[gene==GENE,minCSlevel_top_component:=cumulative_PIP[as.matrix(QTL_perm[gene==GENE,.(1:.N,top_component)])]]
    QTL_perm[gene==GENE,pip_top_component:=susie.perm$alpha[as.matrix(QTL_perm[gene==GENE,.(top_component,1:.N)])]]
    QTL_perm[gene==GENE,lbf_top_component:=susie.perm$lbf[top_component]]
    QTL_perm[gene==GENE,converged:=rep(susie.perm$converged,.N)]
  }else{
    # if <50 variants, set FineMapping values to NA and converged to FALSE
    # for observed data
    QTL_assoc[gene==GENE,rank_top_component:=NA]
    QTL_assoc[gene==GENE,minCSlevel_top_component:=NA]
    QTL_assoc[gene==GENE,pip_top_component:=NA]
    QTL_assoc[gene==GENE,lbf_top_component:=NA]
    QTL_assoc[gene==GENE,converged:=FALSE]
    # for permuted data
    QTL_perm[gene==GENE,rank_top_component:=NA]
    QTL_perm[gene==GENE,minCSlevel_top_component:=NA]
    QTL_perm[gene==GENE,pip_top_component:=NA]
    QTL_perm[gene==GENE,lbf_top_component:=NA]
    QTL_perm[gene==GENE,converged:=FALSE]
  }
  toc()
}
converged=QTL_assoc[,.(converged=any(converged)),by=gene]
print(converged[,paste('FineMapping converged for',sum(converged),'/',length(converged),'genes')])

fwrite(QTL_assoc,file=sprintf('%s/eQTL_FineMapped_chr%s.txt.gz',FINEMAP_DIR,CHR),sep='\t')
fwrite(QTL_perm,file=sprintf('%s/eQTL_FineMapped_chr%s_permuted.txt.gz',FINEMAP_DIR,CHR),sep='\t')
