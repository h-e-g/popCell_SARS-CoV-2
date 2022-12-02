################################################################################
################################################################################
# File name: 2a1__popDEGs_popDRGs.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use linear models to estimate population effects on immune response
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./2a__popDEGs_popDRGs__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))
OUT_DIR=DAT_POPDIFF_DIR

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))
source(sprintf("./misc_plots.R",MISC_DIR))

# define default values
CELLTYPE='lineage'
COV_RUN_NAME='lineage_condition__CellPropLineage_noSVs'
GET_LOGFC=FALSE
PERM=0 # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)

RUN_ID='220409' # date at which the script is launched, used to identify the results unambiguously

STATE='condition'
NLIBS=125 # number of libraries used
ncell_threshold=1 # minimum number of cells to consider CPM in an individual

COV_DIR='1__transcriptome_processing/Covariates/'
EXPRESSION_FILE=NULL

# update parameter values based on arguments provided in 2a1__popDEGs_popDRGs.sh
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use. Will be used for naming of output files
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COV_DIR containing the covariates)
	if (cmd[i]=='--perm' | cmd[i]=='-p' ){PERM = cmd[i+1]} # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)
	if (cmd[i]=='--runid' | cmd[i]=='-d' ){RUN_ID = cmd[i+1]} # ID of the run (e.g., date)
	if (cmd[i]=='--logfc' | cmd[i]=='-g' ){GET_LOGFC = as.logical(cmd[i+1])}
}


EXPRESSION_FILE=sprintf('1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz',NLIBS,CELLTYPE,STATE)

# define run name
RUN_NAME=sprintf('popDiff___%s___%s',COV_RUN_NAME,RUN_ID)


# create output directories
dir.create(sprintf('%s/%s',OUT_DIR,RUN_NAME))
dir.create(sprintf('%s/%s/perm%s',OUT_DIR,RUN_NAME,PERM))

################################################################################
# Transcript-per-million matrix generation and annotation

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]

# remove donors with < 500 cells in at least one condition
filtered_IIDs=fread("1__transcriptome_processing/data/low_cellcount_donors.tsv",sep='\t',header=F)[,V1]
Expression=Expression[!IID%chin%filtered_IIDs,]

# consider only differences between Central Africans and West Europeans
Expression=Expression[POP!="ASH",]
toc()

# if estimating population effects on the transcriptional responses
# compute log-fold changes (logFC) in gene expression relative to NS
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
	Expression[,logFC:=logCPM-logCPM.NS]
}

# define list of celltypes/lineage to consider
groups=Expression[,.N,by=.(celltype,state)]
used_IIDs=unique(Expression[,.(ind=IID,POP)])
#

# randomize individuals
if(PERM>0){
	set.seed(PERM)
	indiv_to_use[,POP:=sample(POP)]
}

# load covariates
CovariatesList=list()
for (i in 1:groups[,.N]){
	myCELLTYPE=groups[i,celltype]
	mySTATE=groups[i,state]
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=fread(sprintf('%s/%s/Covariates__%s_%s.tsv.gz',COV_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE))
}

################################################################################
# Estimate population effects on gene expression/response

# declare function
estimate_population_effects=function(DT, CovariatesList,testID){
  # print test identifier
	cat(testID,'\n')

	# load data
	Covariates=CovariatesList[[gsub('(.*)__(.*)','\\1__\\2',testID)]]
	DT=merge(DT,Covariates, by='IID')
	DT[,POP:=indiv_to_use[match(IID,ind),POP]]

  # define set of covariates to use
	notUsed=c('IID','Symbol','celltype','state','ID','GenderF')
	toUse=c('logCPM','POP',colnames(Covariates[,!c('IID','GenderF')]))

  # build model
	model_ALL=lm(logCPM~ I(POP=='EUB') + . ,data=DT[,..toUse])

	# compute Sandwich estimators (see Methods)
	coefs_ALL=coeftest(model_ALL, vcovHC(model_ALL, type = "HC3")) # See doi: 10.1080/00031305.2000.10474549 for justification

	# extract and format results
	coefs_ALL=data.table(variable=c('intercept','popEUB',colnames(Covariates[,!c('IID','GenderF')])),coefs_ALL[,],df=as.numeric(attr(coefs_ALL,'df')))
	coefs_ALL=melt(coefs_ALL,variable.name='stat',id.vars='variable')
	coefs_ALL[,model:="ALL"]
	coefs_ALL=coefs_ALL[variable!='intercept',]
	coefs_ALL
}

# apply function
tic('peforming popDE analysis')
popDiff_effects=Expression[,
	estimate_population_effects(.SD,CovariatesList,testID=paste(celltype,state,Symbol,sep='__')),
	keyby=.(celltype,state,ID,Symbol)
	]
toc()

# reshape results
tic('spread summary stats across columns')
popDiff_effects_wide=dcast(popDiff_effects,celltype+state+ID+Symbol+variable+model~stat)
setnames(popDiff_effects_wide,c("Estimate","Std. Error","t value","Pr(>|t|)"),c('beta','se','t.value','p.value'))
toc()

# compute FDR for each test separately across all cell types and conditions
tic('compute FDR for each test and write output')
popDiff_effects_wide[se>0,FDR:=p.adjust(p.value,'fdr'),by=.(variable,model)]
toc()

tic('extract population difference parameters')
popDiff=popDiff_effects_wide[ (model=='ALL'&variable=='popEUB', .(ID,Symbol,model,variable,celltype,state,beta,se,p.value,FDR)]
popDiff[,comp:=case_when(model=='ALL'&variable=='popEUB'~'EUB_ref_AFB',TRUE~'other')]
popDiff[is.na(se),se:=0]
toc()

fwrite(popDiff,sprintf("2__population_differences/data/%s_%s.tsv",
	RUN_NAME,	ifelse(PERM>0,sprintf("perm",PERM),"perm0")), sep="\t")
