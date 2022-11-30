################################################################################
################################################################################
# File name: 1d1__compute Covariates.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Computation of covariates (celltype proportions, SVs, etc)
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./1d1__computeCovariates.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))

# read-in parameters
args <- commandArgs(TRUE)
# set default values

# celltype frequency per individual and condition

CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
EXPRESSION_FILE=NULL # expression file to be used.
COV_RUN_NAME=NULL # name of the folder where the Covariates files are saved
SVA=0 # default no SV
GET_LOGFC=FALSE # should this be done from logFC or gene Expression levels

# fixed parameters
NLIBS=125 # number of libraries in the dataset (for naming)
STATE='condition' # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
ncell_threshold=1 # minimal number of cells to consider a sample

CELL_PROP=TRUE # default :include cellular proportions to covariates
CELL_PROP_TYPE=CELLTYPE # resolution at which cellular proportions should be computed

# update parameter values based on provided arguments
for (i in 1:length(args)){
  if (args[i]=='--expression' | args[i]=='-m' ){EXPRESSION_FILE = args[i+1]} # path to expression file used to compute SVs with one line per IID, CELLTYPE, STATE, and gene
  if (args[i]=='--celltype' | args[i]=='-t' ){CELLTYPE = args[i+1]} # celltype variable to use. Will be used for naming of output files
	if (args[i]=='--prop' | args[i]=='-l' ){CELL_PROP = args[i+1]} # should cellular proportions be included (% cells from each celltype/lineage)
  if (args[i]=='--proptype' | args[i]=='-l' ){CELL_PROP_TYPE = args[i+1]} # should cellular proportions included be considered at celltype or lineage level
	if (args[i]=='--sv' | args[i]=='-v' ){SVA =  as.numeric(args[i+1])} # number of Surrogate variables to include (-1: determinaed automatically, 0: none, n>0 : only use the first n SVs )
	if (args[i]=='--run_id' | args[i]=='-r' ){COV_RUN_NAME = args[i+1]} # output directory
	if (args[i]=='--logfc' | args[i]=='-g' ){GET_LOGFC = as.logical(args[i+1])}
}

######## objective of the script:
# for each level of CELLTYPE and STATE, compute relevants covariates to use for popDE/eQTL.
# batch corrected Expression object generated in script 1c2
# Cellcounts folder generated in script 1c1 with cell types proportions for each individual
# define the celltype variable to use when adjusting on cell proportions

if (is.null(EXPRESSION_FILE)){
    EXPRESSION_FILE=sprintf('%s/BatchAdjusted_logCPM_%slibs__per_%s_%s_annotated.tsv.gz',DATA_DIR,NLIBS,CELLTYPE,STATE)
}

# if(is.null(COV_RUN_NAME)){
#         RUN_LOGFC=ifelse(GET_LOGFC,'_logFC','')
# 	RUN_CELL_PROP=ifelse(CELL_PROP==TRUE,'_CellProp','')
#         if(!is.null(CELLTYPE_PROP)){
# 	RUN_CELLTYPE_PROP=sprintf("%s%s",toupper(substr(CELLTYPE_PROP,1,1)),substr(CELLTYPE_PROP,2,nchar(CELLTYPE_PROP)))
#         } else {RUN_CELLTYPE_PROP=""}
# 	RUN_SUBSET=ifelse(SUBSET_PROP==TRUE,'_Subsets','')
# 	RUN_SVs=ifelse(SVA>0,sprintf('_%sSVs',SVA),ifelse(SVA<0,'_SVs',''))
#
# 	COV_RUN_NAME=sprintf('%s_%s%s_%s%s%s%s',CELLTYPE,STATE,RUN_LOGFC,RUN_CELL_PROP,RUN_CELLTYPE_PROP,RUN_SUBSET,RUN_SVs)
# }

############################################
####		  Compute SVs & covariates 		 #####
############################################
MORTALITY_FILE=sprintf('%s/single_cell/project/pop_eQTL/Cell_count_library_mortality.txt',EVO_IMMUNO_POP_ZEUS)
CellMortality=fread(MORTALITY_FILE)
setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
mod=CellMortality[,lm(Mortality~CellCount)]
CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)

if(GET_LOGFC==TRUE){
  Expression_NS=Expression[state=="NS",]
  if(CELLTYPE=="celltype"){
    Expression_NS_CD14.INFECTED=Expression_NS[celltype=='MONO.CD14',]
    Expression_NS_CD14.INFECTED[,celltype:='MONO.CD14.INFECTED']
    Expression_NS=rbind(Expression_NS,Expression_NS_CD14.INFECTED)
    }
	Expression=Expression[state!="NS",]
	# obtain logFC
  Expression=merge(Expression,Expression_NS,by=c('IID','celltype','ID','Symbol','POP'),suffix=c('','.NS'))
	Expression[,logCPM:=logCPM-logCPM.NS]
}
toc()

dir.create(sprintf('%s/Covariates/',OUT_DIR))
dir.create(sprintf('%s/Covariates/%s/',OUT_DIR,COV_RUN_NAME))
dir.create(sprintf('%s/Covariates/%s/SVs/',OUT_DIR,COV_RUN_NAME))

celltype_lineage=fread(sprintf("%s/lineages_celltype.tsv",DATA_DIR))
# ref_celltype="T.CD4.N"

tic('computing SVs and covariates')
# loops across celltypes/conditions
groups=Expression[,.N,keyby=.(celltype,state)]
for (i in 1:groups[,.N]){
  ERROR=FALSE
  myCELLTYPE=groups[i,celltype]
  mySTATE=groups[i,state]
	Expression_select = Expression[celltype==myCELLTYPE & state==mySTATE & ncells >= ncell_threshold,]
#
  Covariates=unique(Expression_select[,.(IID,Age,POP,Gender, ncells)])
  Covariates[,Gender:=relevel(as.factor(Gender),'M')]
  Covariates$Mortality=CellMortality[,setNames(Mortality,IID)][Covariates$IID]

  Cov_mat=model.matrix(~POP+ncells+Gender+Age+Mortality,Covariates)
  rownames(Cov_mat)=Covariates$IID

  if (CELL_PROP==TRUE) {
    if(CELLTYPE_PROP=="lineage"){
      # add percentage of various cell types into the model
      cellProps=fread(sprintf('1__transcriptome_processing/Cellcounts/Pct_celltype_by_IID_%s.tsv.gz',STATE))
      id.vars=c('IID','state',myCELLTYPE)
      #cellTypes_toKeep=setdiff(colnames(cellProps),id.vars)
      if (CELLTYPE=='lineage'){
        lng=myCELLTYPE
        cellTypes_toKeep=celltype_lineage[lineage==lng,celltype]
      } else if (CELLTYPE=='celltype') {
        lng=celltype_lineage[celltype==myCELLTYPE,lineage]
        cellTypes_toKeep=celltype_lineage[lineage==lng,celltype]
      }
      # if (grepl("MONO",myCELLTYPE)|grepl("OTHER",myCELLTYPE)){
      #       cellTypes_toKeep=cellTypes_toKeep[-which(grepl("INFECTED",cellTypes_toKeep))]
      # }
      #     iid_merge=match(rownames(Cov_mat),cellProps$IID)
      #     Cov_mat=cbind(Cov_mat,cellProps[match(paste0(rownames(Cov_mat),mySTATE),paste0(IID,state)),mget(cellTypes_toKeep)])
      # if ((grepl("MONO",myCELLTYPE)|grepl("OTHER",myCELLTYPE))&mySTATE=="IAV"){
      #     rownames(Cov_mat)=cellProps$IID[iid_merge]
      #     Cov_mat=cbind(Cov_mat,cellProps[state=="IAV",][match(rownames(Cov_mat),IID),"MONO.CD14.INFECTED"])
      # }
    } else if (CELLTYPE_PROP=="celltype"){
      cellProps=fread(sprintf('%s/Cellcounts/Pct_celltype_by_IID_%s.tsv.gz',DATA_DIR,STATE))
      id.vars=c('IID','state',myCELLTYPE)
      cellTypes_toKeep=setdiff(colnames(cellProps),id.vars)
      attr(cellTypes_toKeep,"names")=cellTypes_toKeep
      if (CELLTYPE=='lineage'){
        # remove reference cell type to avoid sum of frequencies equal to 1
        cellTypes_toKeep=cellTypes_toKeep[-which(cellTypes_toKeep==ref_celltype)]
      }
    }
    if(mySTATE!='IAV'){
      cellTypes_toKeep=cellTypes_toKeep[-which(grepl("INFECTED",cellTypes_toKeep))]
    }
    Cov_mat=cbind(Cov_mat,cellProps[state==mySTATE,][match(rownames(Cov_mat),IID),mget(cellTypes_toKeep)])
    zeroCols=names(which(apply(Cov_mat==0,2,all)))
    Cov_mat=Cov_mat[,mget(colnames(Cov_mat)[!colnames(Cov_mat)%in%zeroCols])]
    if(all(Cov_mat$GenderF==Cov_mat$POPASH)){

      Cov_mat[,GenderF:=NULL]
    }
  }
  # if(SUBSET_PROP==TRUE){
  #   subsetProps=fread(sprintf('%s/Cellcounts/Pct_cluster_by_%s_%s/Pct_cluster_by_IID_in_%s_%s.tsv.gz',DATA_DIR,CELLTYPE,STATE,myCELLTYPE,mySTATE))
  #   id.vars=c('IID','state',myCELLTYPE)
  #   meanPct=apply(cellProps[state==mySTATE, !..id.vars],2,mean,na.rm=T)
  #   cellTypes_toKeep=names(rev(sort(meanPct))[-1]) # we use the most frequent cell type as reference.
  #   Cov_mat=cbind(Cov_mat,subsetProps[match(rownames(Cov_mat),IID),mget(cellTypes_toKeep)])
  # }
  if(SVA!=0){
  	# prepare TPM matrix for SV computations
  	TPM_DT=dcast(Expression_select, ID + Symbol~IID,value.var="logCPM")
    id.vars=c('ID','Symbol')
  	TPM_Mat=as.matrix(TPM_DT[,!..id.vars])
  	rownames(TPM_Mat)=TPM_DT[,ID]
  	TPM_Mat=TPM_Mat[,Covariates$IID]
    add.noise=function(Matrix,epsilon=1e-7){
        nr=nrow(Matrix)
        nc=ncol(Matrix)
        Matrix = Matrix + matrix(rnorm(nc*nr,0,epsilon),nr,nc)
        Matrix
    }
  	TPM_Mat= add.noise(TPM_Mat,1e-7) # add some noise to avoid numerical issues

  	#compute SVs
    if(SVA>0){
       print(sprintf("Computing SVs for %s_%s:",myCELLTYPE,mySTATE))
  	   SVs=sva(TPM_Mat,mod=as.matrix(Cov_mat),n.sv=SVA,method="two-step")
    }else{
       print(sprintf("Computing SVs for %s_%s:",myCELLTYPE,mySTATE))
         SVs=try(sva(TPM_Mat,mod=as.matrix(Cov_mat),method="two-step"))
    }
    if(class(SVs)!='try-error'){
        rownames(SVs$sv)=Covariates$IID
        if(ncol(SVs$sv)>0){
  	       colnames(SVs$sv)=paste('SV',1:ncol((SVs$sv)), sep='')
           }
        fwrite(data.table(IID=Covariates$IID,SVs$sv),file=sprintf('%s/Covariates/%s/SVs/SVs__%s_%s.tsv.gz',OUT_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE),sep='\t')
  	    if(SVA>0){
          SV_mat=SVs$sv[,1:min(ncol(SVs$sv),SVA)]
        }else{
  		    SV_mat=SVs$sv
  	    }
  	     Cov_mat=cbind(Cov_mat,SV_mat)
       }else{
         ERROR=TRUE
         print(sprintf("Error in SV computation for %s_%s.",myCELLTYPE,mySTATE))
       }
  }
 if(ERROR==FALSE){
  colnames(Cov_mat)=make.names(colnames(Cov_mat))
  COVAR_NAME=setdiff(colnames(Cov_mat),'POPAFB')
  Covariates=data.table(IID=Covariates$IID,Cov_mat[,-1])
  fwrite(Covariates,file=sprintf('%s/Covariates/%s/Covariates__%s_%s.tsv.gz',OUT_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE),sep='\t')
  }
 }
toc()
