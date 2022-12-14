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
EVO_IMMUNO_POP_ZEUS = '/pasteur/zeus/projets/p02/evo_immuno_pop'
OUT_DIR=DAT_POPDIFF_DIR

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))
source(sprintf("./misc_plots.R",MISC_DIR))

# define default values
CELLTYPE='lineage'
STATE='condition'
ADD_CELLPROP=FALSE
PERM=0 # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)
NLIBS=133 # number of libraries used
RUN_ID='210131' # date at which the script is launched, used to identify the results unambiguously
ncell_threshold=1 # minimum number of cells to consider CPM in an individual
COV_DIR=NULL
COV_RUN_NAME='Covariate_CT8_condition_nosubset_noprop_proptype_noSV'
GET_LOGFC=FALSE

META_DATA_FILE="../../0__barcode_processing/data/experimental_metadata.tsv"
MORTALITY_FILE="../../0__barcode_processing/data/per_library_mortality.txt"
EXPRESSION_FILE=NULL

# update parameter values based on arguments provided in 2a1__popDEGs_popDRGs.sh
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use. Will be used for naming of output files
	if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
	if (cmd[i]=='--covdir' | cmd[i]=='-r' ){COV_DIR = cmd[i+1]} # directory where covariates can be found (default: single_cell/project/pop_eQTL/Covariates )
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COV_DIR containing the covariates)
	if (cmd[i]=='--perm' | cmd[i]=='-p' ){PERM = cmd[i+1]} # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output dir : folder to save popDiff result tables (default single_cell/project/pop_eQTL/popDE )
	if (cmd[i]=='--runid' | cmd[i]=='-d' ){RUN_ID = cmd[i+1]} # ID of the run (e.g., date)
	if (cmd[i]=='--nlibs' | cmd[i]=='-n' ){NLIBS = cmd[i+1]} # number of libraries used
	if (cmd[i]=='--cellprop' | cmd[i]=='-c' ){ADD_CELLPROP = as.logical(cmd[i+1])} # cellprop : should cell proportions (computed as mean of the 3 conditions) be included in the model: default : FALSE
	if (cmd[i]=='--figdir' | cmd[i]=='-f' ){FIGURE_DIR = cmd[i+1]} # figure dir : where ouput figures should be saved (default single_cell/project/pop_eQTL/figurespopDE )
	if (cmd[i]=='--expression' | cmd[i]=='-e' ){EXPRESSION_FILE = cmd[i+1]}
	if (cmd[i]=='--metadata' | cmd[i]=='-m' ){META_DATA_FILE = cmd[i+1]}
	if (cmd[i]=='--logfc' | cmd[i]=='-g' ){GET_LOGFC = as.logical(cmd[i+1])}
}

DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences",EVO_IMMUNO_POP_ZEUS)
if(is.null(EXPRESSION_FILE)){
	EXPRESSION_FILE=sprintf('%s/BatchAdjusted_logCPM_%slibs__per_%s_%s_annotated.tsv.gz',DATA_DIR,NLIBS,CELLTYPE,STATE)
}

if(is.null(COV_DIR)){
	COV_DIR= sprintf("%s/Covariates",DATA_DIR)
}

# define run name
RUN_LOGFC=ifelse(GET_LOGFC,'_logFC','')
RUN_CELLPROP=ifelse(ADD_CELLPROP,'_19cellProp','')
RUN_NAME=sprintf('%s_%s_AFBEUB_%s_%s%s_%s',CELLTYPE,STATE,RUN_LOGFC,COV_RUN_NAME,RUN_CELLPROP,RUN_ID)


# create output directories
dir.create(sprintf('%s/%s',OUT_DIR,RUN_NAME))
dir.create(sprintf('%s/%s/perm%s',OUT_DIR,RUN_NAME,PERM))
dir.create(sprintf('%s/popDE/%s',FIGURE_DIR,RUN_NAME))
dir.create(sprintf('%s/popDE/%s/perm%s',FIGURE_DIR,RUN_NAME,PERM))

# load sample meta data
CellMortality=fread(MORTALITY_FILE)
setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
mod=CellMortality[,lm(Mortality~CellCount)]
CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

meta_data=fread(META_DATA_FILE, sep='\t')

################################################################################
# Transcript-per-million matrix generation and annotation

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]

# remove donors with < 500 cells in at least one condition
filtered_IIDs=fread("../../1__transcriptome_processing/data/low_cellcount_donors.tsv",sep='\t',header=F)[,V1]
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

# define list of clusters to consider
groups=Expression[,.N,by=.(celltype,state)]
used_IIDs=unique(Expression[,.(ind=IID,POP)])

SELECTED_CELLTYPES=unique(groups[,celltype])
if(CELLTYPE=='lineage'){
  SELECTED_CELLTYPES=SELECTED_CELLTYPES[!grepl("^OTHER",SELECTED_CELLTYPES)]
  groups=groups[celltype!='OTHER',]
}

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
	# add cell mortality
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]][,-c('POPASH','POPEUB','ncells','CellCount')]
	# add aggregated cell proportions
	if(ADD_CELLPROP){
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=merge(CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]],CellPct_mean,by='IID')
	}
}

################################################################################
# Estimate population effects on gene expression/response

# declare function
estimate_population_effects=function(DT, CovariatesList,testID){
  # print test identifier
	cat(testID,'\n')

	# load data
	Covariates=CovariatesList[[gsub('(.*)__(.*)__(.*)','\\1__\\2',testID)]]
	DT=merge(DT,Covariates, by='IID')
	DT[,POP:=indiv_to_use[match(IID,ind),POP]]

  # define set of covariates to use
	notUsed=c('IID','Symbol','celltype','state','ID','GenderF')
	toUse=c('logCPM','POP',colnames(Covariates[,!c('IID','GenderF')]))

  # build model
	model_ALL=lm(logCPM~ I(POP=='EUB') + . ,data=DT[,..toUse]) # +I((POP=='AFB') - (POP=='EUB')

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
popDiff_effects=Expression[
  celltype%chin%SELECTED_CELLTYPES,
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

fwrite(popDiff,sprintf("2__population_differences/data/popdiff__AFB_EUB__%slineage_condition%s__%s%s_%s.tsv",
  ifelse(GET_LOGFC,"logFC_",""),
	ifelse(GET_LOGFC,"_logFC",""),
	ifelse(ADD_CELLPROP,"CellPropLineage",""),
	RUN_ID,
	ifelse(PERM>0,sprintf("perm",PERM),"perm0")),
sep="\t")

################################################################################
# Assess significance of effect size decrease after adjustments

# which population differences are actually lost (i.e., associated to signif.
# decrease in effect size) after adjusting for cellular composition?

# no adjustment > within-lineage (lineage) adjustment
# x := larger group (i.e., no)
# y := smaller group (i.e., lineage)

# estimate correlation between raw (beta_x) and adjusted (beta_y) under the null (permuted data)
# even if the data are permuted, the beta_x and beta_y are estimated from the same data and could hence be correlated
# this dependence needs to be accounted for when performing the t-test, lest we underestimate the number of significant changes

# load raw population effect estimates on gene expression and response
ANALYSE="popdiff__AFB_EUB__lineage_condition__220409_perm1.tsv"
popDE_noCellProp_null=fread(sprintf("2__population_differences/%s",ANALYSE))
popDE_noCellProp_null[,type:="expression"]
ANALYSE="popdiff__AFB_EUB__logFC_lineage_condition_logFC__220409_perm1.tsv"
popDR_noCellProp_null=fread(sprintf("2__population_differences/%s",ANALYSE))
popDR_noCellProp_null[,type:="response"]

# load adjusted population effect estimates on gene expression and response
ANALYSE="popdiff__AFB_EUB__lineage_condition__CellPropLineage_220409_perm1.tsv"
popDE_lineageCellProp_null=fread(sprintf("2__population_differences/%s",ANALYSE))
popDE_lineageCellProp_null[,type:="expression"]
ANALYSE="popdiff__AFB_EUB__logFC_lineage_condition_logFC__CellPropLineage_220409_perm1.tsv"
popDR_lineageCellProp_null=fread(sprintf("2__population_differences/%s",ANALYSE))
popDR_lineageCellProp_null[,type:="response"]

# annotate, merge and reshape
popDE_noCellProp_null[,cellprop_adjust:="no"]
popDE_lineageCellProp_null[,cellprop_adjust:="lineage"]
popDE_adjusted_null=rbind(popDE_noCellProp,popDE_lineageCellProp)
popDE_adjusted_null[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]
popDE_adjusted_null_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

popDR_noCellProp_null[,cellprop_adjust:="no"]
popDR_lineageCellProp_null[,cellprop_adjust:="lineage"]
popDR_adjusted_null=rbind(popDR_noCellProp,popDR_lineageCellProp)
popDR_adjusted_null[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]
popDR_adjusted_null_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

popDE_adjusted_null_wide=dcast(popDE_adjusted_null_comp,idx~group,value.var=c("beta","FDR","se"))

# correlation between raw and adjusted population effects on expression under null
corr_beta_no_lineage_de=popDE_adjusted_null_wide[,cor.test(beta_x,beta_y)]$estimate
# x := no, y := lineage
# > popDE_adjusted_null_wide[,cor.test(beta_x,beta_y)]
#
#         Pearson's product-moment correlation
#
# data:  beta_x and beta_y
# t = 818.66, df = 190063, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8816486 0.8836351
# sample estimates:
#       cor
# 0.8826458

# test for significant differences in raw and adjusted effect sizes, given
# correlation under the null
significant_popde_losses=lapply(c("DE","DR"),function(exp_resp){
  noCellProp=get(sprintf("pop%s_noCellProp",exp_resp))
  lineageCellProp=get(sprintf("pop%s_lineageCellProp",exp_resp))
  noCellProp[,cellprop_adjust:="no"]
  lineageCellProp[,cellprop_adjust:="lineage"]
  adjusted=rbind(noCellProp,lineageCellProp)
  adjusted[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]

  adjusted_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

  signif_in_xy=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[(FDR_x<0.01&abs(beta_x)>0.2)&(FDR_y<0.01&abs(beta_y)>0.2),idx]
  signif_in_x=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[(FDR_x<0.01&abs(beta_x)>0.2)&!(FDR_y<0.01&abs(beta_y)>0.2),idx]
  signif_in_y=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[!(FDR_x<0.01&abs(beta_x)>0.2)&!(FDR_y<0.01&abs(beta_y)>0.2),idx]

  adjusted_comp$signif_in=case_when(
    adjusted_comp$idx%in%signif_in_xy~"xy",
    adjusted_comp$idx%in%signif_in_x~"x",
    adjusted_comp$idx%in%signif_in_y~"y"
  )

  adjusted_wide=dcast(adjusted_comp[!is.na(signif_in)],idx+signif_in~group,value.var=c("beta","FDR","se"))

  adjusted_wide=adjusted_wide[!grepl("^IAV|^SARS",idx),]


  betas_t=adjusted_wide[,.(t=(beta_x-beta_y)/(sqrt(se_x^2+se_y^2-2*corr_beta_no_lineage_de*se_x*se_y))),by=.(idx)]

  adjusted_wide$t=betas_t[,setNames(t,idx)][adjusted_wide$idx]

  adjusted_wide$different=ifelse(abs(adjusted_wide$t)>1.96,"Different","Same")

  adjusted_wide$Symbol=adjusted[,setNames(Symbol,idx)][adjusted_wide$idx]
  adjusted_wide$state=str_split(adjusted_wide$idx,"_",simplify=T)[,3]
  adjusted_wide$celltype=str_split(adjusted_wide$idx,"_",simplify=T)[,2]

  adjusted_labels=adjusted_wide[different=="Different"&signif_in!="all",][order(-abs(t)),][,head(.SD,3),by=.(signif_in,state,celltype)][,.(beta_x,beta_y,signif_in,celltype,state,Symbol)]

  adjusted_wide[,difference_type:=case_when(signif_in=="x"~"difference_lost",signif_in=="y"~"difference_gained",signif_in=="xy"~"difference_kept")]
  setnames(adjusted_wide,"different","difference_signif")
  result=adjusted_wide[,-"signif_in"]
	result[,type:=exp_resp]
  return(result)
})%>%rbindlist()

fwrite(significant_popde_losses,"2__population_differences/data/cell_composition_adjustment_beta_changes.tsv",sep="\t")
