######## objective of the script:
# for each level of CELLTYPE and STATE, compute population differences  on observed/permuted data

######## output of the script:
# a folder will be created in OUT_DIR with a name that reflects the set of covariates used (COV_RUN_NAME), whether analyses are adjusted on sex and if results where combined across populations
# in this folder (call RUN_NAME thereafter), one sub folder per permutation wil be done.
# PopDiffResults__per_{CELLTYPE}_{STATE}.tsv.gz contains all results from the two population difference test models (including covariates)

# model : ALL Expr ~ popASH + popEUB + Covariates
# model : DIF Expr ~ popAFB + popASH + Covariates

# celltype	state   ID                Symbol  variable model   beta    se        t.value   p.value   df    FDR
# Bcells    ACTIVE  ENSG00000000419   DPM1    Age      ALL   -0.0019   0.0025    -0.7746   0.4395   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    Age      DIF   -0.0019   0.0025    -0.7746   0.4395   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    Gender   ALL   -0.0099   0.0993    -0.1006   0.9199   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    Gender   DIF   -0.0099   0.0993    -0.1006   0.9199   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    ncells   ALL    0.0001   0.0001     0.6241   0.5333   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    ncells   DIF    0.0001   0.0001     0.6241   0.5333   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    popAFB   DIF   -0.0354   0.0655    -0.5411   0.5890   191    0.9847
# Bcells    ACTIVE  ENSG00000000419   DPM1    popASH   ALL    0.0485   0.0893     0.5432   0.5875   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    popASH   DIF    0.0130   0.0792     0.1648   0.8692   191    0.9999
# Bcells    ACTIVE  ENSG00000000419   DPM1    popEUB   ALL    0.0354   0.0655     0.5411   0.5890   191    0.9847
# Bcells    ACTIVE  ENSG00000000457   SCYL3   Age      ALL   -0.0004   0.0029    -0.1707   0.8646   191    0.9999
# Bcells    ACTIVE  ENSG00000000457   SCYL3   Age      DIF   -0.0004   0.0029    -0.1707   0.8646   191    0.9999
# Bcells    ACTIVE  ENSG00000000457   SCYL3   Gender   ALL   -0.0284   0.0926    -0.3073   0.7589   191    0.9999

# mash_fitted_${CELLTYPE}_${STATE}_3pops.RDS contain the fitted mash model object

# popDiffResults_with_mash.tsv.gz contains all population difference results for EUBvsAFB, ASHvsAFB & ASHvsEUB, with both raw results and mashr estimated parameteres (beta_mash, se_mash and lfsr_mash)
# ID                Symbol  model   variable    celltype       state    beta      se       p.value   FDR       comp         beta_mashr    se_mashr    lfsr_mashr
# ENSG00000000419   DPM1    ALL    popASH       B.M         ACTIVE     0.1040    0.0984    0.2919    0.999    ASH_ref_AFB   -1.3e-05      0.0007      0.9992
# ENSG00000000419   DPM1    ALL    popASH       B.M         RESTING    0.4665    0.2438    0.0573    0.999    ASH_ref_AFB   -1.5e-05      0.0008      0.9992
# ENSG00000000419   DPM1    ALL    popASH       B.M         SHARED    -0.1414    0.1933    0.4660    0.999    ASH_ref_AFB   -1.1e-05      0.0006      0.9992
# ENSG00000000419   DPM1    ALL    popASH       B.N         ACTIVE    -0.0492    0.1334    0.7129    0.999    ASH_ref_AFB   -1.4e-05      0.0008      0.9992

# Nb_popDE.tsv contains number of popDE for various thresholds.

# some figures (Nb of popDE, Sharing between pops, sharing between conds) are also outputted in {FIGURE_DIR}/{RUN_NAME}

###### requirements:
### define CELLTYPE and STATE
# ideally CELLTYPE and CELLSTATE, should not contain any '_'
# eg. CELLTYPE='celltype.19level'
#     STATE='activation.state'
#
# or CELLTYPE='celltype.8level'
#     STATE='condition'

### batch corrected Expression object generated in script 6c

### Covariates folder generated in script 6e

### set NLIBS to the number of libraries in the current iteration (for naming purposes)


##################### Expected inputs:
# Expression file: with logCPM of all genes, per IID x celltype x state
# metadata file: with 1 line per IID x COND x RUN, used to annotate samples

# define default values
CELLTYPE='lineage'
STATE='condition'
ADD_CELLPROP=FALSE
PERM=0 # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)
NLIBS=125 # number of libraries used
RUN_ID='210131' # date at which the script is launched (ar any character string), is used to identify the results unambiguously.
ncell_threshold=1 # minimum number of cells to consider CPM in an individual
COV_DIR=NULL
COV_RUN_NAME='Covariate_CT8_condition_nosubset_noprop_proptype_noSV'
GET_LOGFC=FALSE

EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
FIGURE_DIR = sprintf("%s/single_cell/project/pop_eQTL/figures/2_population_differences",EVO_IMMUNO_POP_ZEUS)
OUT_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences/popDE",EVO_IMMUNO_POP_ZEUS)

META_DATA_FILE=sprintf('%s/popCell_data/00_CRF/scrnaseq_popbased_metadata_full_long.tsv',EVO_IMMUNO_POP_ZEUS)
MORTALITY_FILE=sprintf('%s/single_cell/project/pop_eQTL/Cell_count_library_mortality.txt',EVO_IMMUNO_POP_ZEUS)
EXPRESSION_FILE=NULL
# update parameter values based on provided arguments
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use. Will be used for naming of output files
	if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
	if (cmd[i]=='--covdir' | cmd[i]=='-r' ){COV_DIR = cmd[i+1]} # directory where covariates can be found (default: single_cell/project/pop_eQTL/Covariates )
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COV_DIR containing the covariates)
	if (cmd[i]=='--perm' | cmd[i]=='-p' ){PERM = cmd[i+1]} # should we permute the population labels (0: no ; n>0: yes, permutation is done with seed n)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output dir : folder to save popDiff result tables (default single_cell/project/pop_eQTL/popDE )
	if (cmd[i]=='--runid' | cmd[i]=='-d' ){RUN_ID = cmd[i+1]} # ID of the run (eg. date)
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

###############@ START SCRIPT
options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggrastr))
suppressMessages(library(ggcorrplot))
suppressMessages(library(scales))
suppressMessages(library(data.table))
# suppressMessages(library(SingleCellExperiment))
# suppressMessages(library(scran))
# suppressMessages(library(dynamicTreeCut))
# suppressMessages(library(BiocNeighbors))
# suppressMessages(library(cowplot))
suppressMessages(library(sandwich))
suppressMessages(library(lmtest))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(sva))

####################
# Define colors

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))
BROAD=ifelse(CELLTYPE=='celltype.8level',TRUE,FALSE)
if(BROAD==TRUE){
	color_cellTypes=color_cellTypes_8level
}else{
	color_cellTypes=color_cellTypes_19level
}
####################
# Define RUN NAME

RUN_LOGFC=ifelse(GET_LOGFC,'_logFC','')
RUN_CELLPROP=ifelse(ADD_CELLPROP,'_19cellProp','')
RUN_NAME=sprintf('%s_%s_AFBEUB_%s_%s%s_%s',CELLTYPE,STATE,RUN_LOGFC,COV_RUN_NAME,RUN_CELLPROP,RUN_ID)

###################
# create output directories
dir.create(sprintf('%s/%s',OUT_DIR,RUN_NAME))
dir.create(sprintf('%s/%s/perm%s',OUT_DIR,RUN_NAME,PERM))

dir.create(sprintf('%s/popDE/%s',FIGURE_DIR,RUN_NAME))
dir.create(sprintf('%s/popDE/%s/perm%s',FIGURE_DIR,RUN_NAME,PERM))

# dir.create(sprintf('%s/%s/Covariates',WORKDIR,RUN_NAME))


##### define sample meta data
CellMortality=fread(MORTALITY_FILE)
setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
mod=CellMortality[,lm(Mortality~CellCount)]
CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

# TODO make sure metadata file is up to date
meta_data=fread(META_DATA_FILE, sep='\t')

###########################################################################################
######### STEP 1: annotate data (Age, Sex, ncells, POP) and generate TPM_Matrix  ##########
###########################################################################################

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]
# remove inds with < 500 cells
MinCell_perCOND=500
filtIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/filteredIID_lessThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
Expression=Expression[!IID%chin%filtIID,]
toc()
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
indiv_to_use=unique(Expression[,.(ind=IID,POP)])

SELECTED_CELLTYPES=unique(groups[,celltype])
if(CELLTYPE=='lineage'){
SELECTED_CELLTYPES=SELECTED_CELLTYPES[!grepl("^OTHER",SELECTED_CELLTYPES)]
groups=groups[celltype!='OTHER',]
}

indiv_to_use=indiv_to_use[POP!='ASH',]

if(PERM>0){
# randomize individuals
	set.seed(PERM)
	indiv_to_use[,POP:=sample(POP)]
	}

# if(ADD_CELLPROP){
# 	meta=fread(sprintf('%s/single_cell/project/pop_eQTL/data/2_population_differences/sce_clean_metadata.tsv',EVO_IMMUNO_POP_ZEUS))
# 	CellPct=meta[,.N,by=.(IID,celltype.19level,COND)]
# 	CellPct=CellPct[,Pct:=N/sum(N),by=.(IID,COND)]
# 	CellPct_mean=CellPct[,.(Pct=mean(Pct)),by=.(IID,celltype.19level)]
# 	CellPct_mean=dcast(CellPct_mean,IID~celltype.19level,fill=0)
# }

CovariatesList=list()
for (i in 1:groups[,.N]){
	myCELLTYPE=groups[i,celltype]
	mySTATE=groups[i,state]
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=fread(sprintf('%s/%s/Covariates__%s_%s.tsv.gz',COV_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE))
	# Add Cell Mortality
#	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=merge(CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]],CellMortality,by='IID')
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]][,-c('POPASH','POPEUB','ncells','CellCount')]
	# Add AGGREGATE Cell PROPORTIONS
	if(ADD_CELLPROP){
	CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]]=merge(CovariatesList[[paste(myCELLTYPE,mySTATE,sep='__')]],CellPct_mean,by='IID')
	}
}

############################################
#### STEP 3:  estimate popDiff Effects #####
############################################
test_POPs=function(DT, CovariatesList,testID){
		cat(testID,'\n')
		Covariates=CovariatesList[[gsub('(.*)__(.*)__(.*)','\\1__\\2',testID)]]
		DT=merge(DT,Covariates, by='IID')
		DT[,POP:=indiv_to_use[match(IID,ind),POP]]

		notUsed=c('IID','Symbol','celltype','state','ID','GenderF')
		toUse=c('logCPM','POP',colnames(Covariates[,!c('IID','GenderF')]))

#		model_ALL=lm(value~I(POP=='ASH') + I(POP=='EUB') + Age + ncells + . ,data=DT[,!..notUsed]) # +I((POP=='AFB') - (POP=='EUB')
		model_ALL=lm(logCPM~ I(POP=='EUB') + . ,data=DT[,..toUse]) # +I((POP=='AFB') - (POP=='EUB')
		coefs_ALL=coeftest(model_ALL, vcovHC(model_ALL, type = "HC3")) # See doi: 10.1080/00031305.2000.10474549 for justification
	  coefs_ALL=data.table(variable=c('intercept','popEUB',colnames(Covariates[,!c('IID','GenderF')])),coefs_ALL[,],df=as.numeric(attr(coefs_ALL,'df')))
		coefs_ALL=melt(coefs_ALL,variable.name='stat',id.vars='variable')
	  coefs_ALL[,model:='ALL']
	  # coefs_ALL=coefs_ALL[variable!='intercept' & !grepl('^SV[0-9]+$',variable),]
		coefs_ALL=coefs_ALL[variable!='intercept',]
		coefs_ALL
	}


tic('peforming popDE analysis')
popDiff_effects=Expression[POP!='ASH' & celltype%in%SELECTED_CELLTYPES,test_POPs(.SD,CovariatesList,testID=paste(celltype,state,Symbol,sep='__')),keyby=.(celltype,state,ID,Symbol)]
# Expression_1=Expression[Symbol %in% Expression[1:300,Symbol],]
# popDiff_effects=Expression_1[POP!='ASH' & celltype%in%SELECTED_CELLTYPES,test_POPs(.SD,CovariatesList,testID=paste(celltype,state,Symbol,sep='__')),keyby=.(celltype,state,ID,Symbol)]
toc()

# fwrite(popDiff_effects,file=sprintf('%s/%s/perm%s/PopDiffResults__per_%s_%s.tsv.gz',OUT_DIR,RUN_NAME,PERM,CELLTYPE,STATE),sep='\t')
#popDiff_effects=fread(sprintf('%s/%s/perm%s/PopDiffResults__per_%s_%s.tsv.gz',OUT_DIR,RUN_NAME,PERM,CELLTYPE,STATE),sep='\t')

############################################
#### STEP 4:  reshape popDiff results  #####
###########################################

tic('spread summary stats across columns')
popDiff_effects_wide=dcast(popDiff_effects,celltype+state+ID+Symbol+variable+model~stat)
setnames(popDiff_effects_wide,c("Estimate","Std. Error","t value","Pr(>|t|)"),c('beta','se','t.value','p.value'))
toc()

tic('compute FDR for each test and write output')
# compute FDR for each test (eg. Age, Gender, PopEUB vs others,  PopAFB vs others, PopEUB vs popAFB) separately across all cell types & activation states)
popDiff_effects_wide[se>0,FDR:=p.adjust(p.value,'fdr'),by=.(variable,model)]
fwrite(popDiff_effects_wide,file=sprintf('%s/%s/perm%s/PopDiffResults__per_%s_%s.tsv.gz',OUT_DIR,RUN_NAME,PERM,CELLTYPE,STATE),sep='\t')
toc()
# popDiff_effects_wide=fread(sprintf('%s/%s/perm%s/PopDiffResults__per_%s_%s.tsv.gz',OUT_DIR,RUN_NAME,PERM,CELLTYPE,STATE))

tic('extract population difference parameters')
###### extract the 3 relevant comparisons
popDiff=popDiff_effects_wide[ (model=='ALL' & (variable=='popEUB' |  variable=='popASH')) | (model=='DIF' & variable=='popASH'), .(ID,Symbol,model,variable,celltype,state,beta,se,p.value,FDR)]
# rename them
popDiff[,comp:=case_when(model=='ALL' & variable=='popEUB'~ 'EUB_ref_AFB',
                        model=='ALL' & variable=='popASH'~ 'ASH_ref_AFB',
                        model=='DIF' & variable=='popASH'~ 'ASH_ref_EUB',
                        TRUE~'other')]
popDiff[is.na(se),se:=0]

print(popDiff[,.N,by=comp])
toc()


########################################################################
####################### 		SETTING UP MASHR 			######################
########################################################################

tic('prepare mash input and estimate null correlations')
	id.vars=c('comp','ID','Symbol')
	group.vars=c('celltype','state')
	beta_mat=dcast(popDiff, comp+ID+Symbol~celltype+state, value.var='beta',fill=0)
	se_mat=dcast(popDiff, comp+ID+Symbol~celltype+state, value.var='se',fill=0)

# generate row and columns ids for beta_mat and se_mat
# (dependent on whether we share info across population comparisons)
cat('where ')
row_IDs=apply(beta_mat[,..id.vars],1,paste,collapse='_')
cat('does ')
prov=apply(popDiff[,..id.vars],1,paste,collapse='_')
popDiff[,row_ID:=prov]
cat('it ')
prov=apply(popDiff[,..group.vars],1,paste,collapse='_')
popDiff[,group:=prov]
cat('crash?')

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

library(mashr)
data.noised = mash_set_data(add.noise(beta_mat[,!..id.vars], Smin, row_IDs),
                              as.matrix(se_mat[,!..id.vars])+Smin)

########################################################################
##################### 				RUNNING MASHR 			######################
########################################################################

# estimate correlations between null tests
Vhat = estimate_null_correlation_simple(data.noised,z_thresh=3)
data.Vhat = mash_update_data(data.noised, V=Vhat)
toc()

tic('estimate correlations between strong signals to define covariance structure between tests under H1')
m.1by1 = mash_1by1(data.Vhat)
strong = get_significant_results(m.1by1, 0.05/ncol(m.1by1$result$lfsr))

if(length(strong)<50){
	WARNING='WARNING: too few strong association to learn covariance structure from, decreasing significance level'
	strong = get_significant_results(m.1by1, 0.05)
	write(WARNING,file=sprintf('%s/%s/perm%s/WARNINGs.txt',OUT_DIR,RUN_NAME,PERM),append=F)
	if(length(strong)<50){
		WARNING='WARNING: still too few strong associations, picking 100 random test'
		write(WARNING,file=sprintf('%s/%s/perm%s/WARNINGs.txt',OUT_DIR,RUN_NAME,PERM),append=T)
		strong=sample(1:nrow(beta_mat),100)
	}
}
data.strong = mash_set_data(add.noise(beta_mat[strong,!..id.vars], Smin, row_IDs[strong]),
                              as.matrix(se_mat[strong,!..id.vars])+Smin, V=Vhat)

# define covariance structure between tests under H1
nCOND=ncol(beta_mat)-length(id.vars)
U.pca = cov_pca(data.strong,min(5,nCOND))
print(names(U.pca))
U.c = cov_canonical(data.Vhat)
toc()

# this crashes, likely due to genes with 0 values for either beta or se.
#U.ed = try(cov_ed(data.strong, U.pca))
#m = mash(data.reduced, Ulist = c(U.ed,U.c), outputlevel = 1)

tic('running mash')
m = mash(data.Vhat, Ulist = c(U.pca,U.c), outputlevel = 1)
m2 = mash(data.Vhat, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2,file=sprintf('%s/%s/perm%s/mash_fitted_%s_%s_AFEUonly.RDS',OUT_DIR,RUN_NAME,PERM,CELLTYPE,STATE))
toc()

tic('get mash results')
lfsr_mashr=get_lfsr(m2)
beta_mashr=get_pm(m2)
se_mashr=get_psd(m2)
toc()

########################################################################
###############			ADDING MASHR RESULTS TO POPDIFF 		################
########################################################################


tic('reshape mash results')
lfsr_mashr_long=melt(data.table(row_ID=row_IDs,lfsr_mashr), id.vars='row_ID', value.name='lfsr_mashr', variable.name='group')
beta_mashr_long=melt(data.table(row_ID=row_IDs,beta_mashr), id.vars='row_ID', value.name='beta_mashr', variable.name='group')
se_mashr_long=melt(data.table(row_ID=row_IDs,se_mashr), id.vars='row_ID', value.name='se_mashr', variable.name='group')

mashr_long=merge(beta_mashr_long,se_mashr_long,by=c('row_ID','group'))
mashr_long=merge(mashr_long,lfsr_mashr_long,by=c('row_ID','group'))

popDiff=merge(popDiff,mashr_long,by=c('row_ID','group'))
popDiff[,group:=NULL]
popDiff[,row_ID:=NULL]
fwrite(popDiff,file=sprintf("%s/%s/perm%s/popDiffResults_with_mash.tsv.gz",OUT_DIR,RUN_NAME,PERM),sep='\t')
toc()

########################################################################
###############  summary plots 		##########################
########################################################################


CountDE=popDiff[,.(FDR1=length(unique(ID[FDR<0.01])),
FDR5=length(unique(ID[FDR<0.05])),
LFSR5=length(unique(ID[lfsr_mashr<0.05])),
LFSR20=length(unique(ID[lfsr_mashr<0.2])),
LFSR5_beta20=length(unique(ID[lfsr_mashr<0.05 & abs(beta_mashr)>0.2])),
LFSR5_beta50=length(unique(ID[lfsr_mashr<0.05 & abs(beta_mashr)>0.5])))
, by=.(celltype,state,comp)]

# ncells=meta[,.N,by=.(IID,celltype.19level,COND)][,.(ncells=mean(ncells)),keyby=.(celltype,state)]
# CountDE=merge(CountDE,ncells,by=c('celltype','state'))
#DT=melt(CountDE,id.vars=c('celltype','state','ncells','comp'))
DT=melt(CountDE,id.vars=c('celltype','state','comp'))

fwrite(DT,file=sprintf("%s/%s/perm%s/Nb_popDE.tsv",OUT_DIR,RUN_NAME,PERM),sep='\t')
DT[,celltype:=make.names(celltype)]

pdf(sprintf("%s/popDE/%s/perm%s/Nb_popDE_freescale.pdf",FIGURE_DIR,RUN_NAME,PERM),height=12)
p <- ggplot(DT,aes(x=paste(celltype,state,sep='-'), y=value, fill=celltype, color=state)) + geom_bar(stat='identity',position='dodge')
p <- p + theme_mary() + facet_grid(variable~comp, scales="free_y")+ xlab('') + ylab('# of popDE genes')
p <- p + scale_fill_manual(values=color_cellTypes) + scale_color_manual(values=color_activation)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

pdf(sprintf("%s/popDE/%s/perm%s/Nb_popDE_sqrt.pdf",FIGURE_DIR,RUN_NAME,PERM),height=12)
p <- ggplot(DT,aes(x=paste(celltype,state,sep='-'), y=value, fill=celltype, color=state)) + geom_bar(stat='identity',position='dodge')
p <- p + theme_mary() + facet_grid(variable~comp)+ xlab('') + ylab('# of popDE genes') + ylim(c(0,10000)) + scale_y_sqrt()
p <- p + scale_fill_manual(values=color_cellTypes) + scale_color_manual(values=color_activation)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

if(PERM==0){
	pi_shared=get_pairwise_sharing(m2)
	DT=as.data.table(expand.grid(1:nrow(pi_shared),1:ncol(pi_shared)))
	colnames(DT)=c('Var1','Var2')
	DT[,sharing:=pi_shared[cbind(Var1,Var2)]]
	DT[,group_1:=colnames(pi_shared)[Var1]]
	DT[,group_2:=colnames(pi_shared)[Var2]]

	# if(COMBINE_POPS){
	# 	DT[,group_1:=gsub('(.*)_(.*)_(ASH|AFB|EUB)_ref_(ASH|AFB|EUB)','\\3_ref_\\4_\\1_\\2',group_1)]
	# 	DT[,group_2:=gsub('(.*)_(.*)_(ASH|AFB|EUB)_ref_(ASH|AFB|EUB)','\\3_ref_\\4_\\1_\\2',group_2)]
	# }

	pdf(sprintf("%s/popDE/%s/Sharing_popDE.pdf",FIGURE_DIR,RUN_NAME),height=10,width=10)
	p <- ggplot(DT,aes(x=group_1, y=group_2, fill=sharing)) + geom_tile()
	p <- p + theme_mary() + scale_fill_distiller(palette="RdYlBu",direction=-1,limits=c(0, 1))
	p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab('')+ylab('')
	print(p)
	dev.off()
}

q('no')
