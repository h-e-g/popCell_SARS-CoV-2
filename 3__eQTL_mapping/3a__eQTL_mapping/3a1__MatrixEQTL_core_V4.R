
##################### Expected inputs:
# Expression file: with logCPM of all genes, per IID x celltype x state
# metadata file: with 1 line per IID x COND x RUN, used to annotate samples

# define default values
TEST=FALSE # if true, run on a samll susbet of SNPs to test code
CELLTYPE='lineage'
STATE='condition'
ncell_threshold=1 # minimum number of cells to consider CPM in an individual
COV_DIR=NULL
ADD_CELLPROP=FALSE

EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
PROJECT = "pop_eQTL"
FIGURE_DIR = sprintf("%s/single_cell/project/%s/figures/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS,PROJECT)
DATA_DIR = sprintf("%s/single_cell/project/%s/data",EVO_IMMUNO_POP_ZEUS,PROJECT)
SCRATCH="/pasteur/appa/scratch/users/mrotival/"
OUT_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
COV_RUN_NAME='lineage_condition__CellPropLineage_SVs'

META_DATA_FILE=sprintf('%s/popCell_data/00_CRF/scrnaseq_popbased_metadata_full_long.tsv',EVO_IMMUNO_POP_ZEUS)
MORTALITY_FILE=sprintf('%s/single_cell/project/pop_eQTL/Cell_count_library_mortality.txt',EVO_IMMUNO_POP_ZEUS)

EXPRESSION_FILE=NULL
RUN_ID="test_run" # <<<<update when rerunning to generate files in a separate folder >>>>

CIS_DIST=1e6 # Distance in Cis to consider in the analysis.
pvCis=1 # min P-value in cis to report a variant
minFreq=0.05 # min allele frequency (per population) to consider a variant
NLIBS=125 # number of libraries considered to the run
GET_LOGFC=FALSE

cmd=commandArgs()
for (i in 1:length(cmd)){
	if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use to define pseudoBulk. Will be used for naming of output files
	if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use to define pseudoBulk (cellular activation state or experimental condition). Will be used for naming of output files
	if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]}
	if (cmd[i]=='--mystate' | cmd[i]=='-o' ){mySTATE = cmd[i+1]} # condition/activationn state for which the eQTL mapping will be done (should be one level of STATE variable)
	if (cmd[i]=='--mycelltype' | cmd[i]=='-y' ){myCELLTYPE = cmd[i+1]} # celltype for which the eQTL mapping will be done (should be one level of CELLTYPE variable)
	if (cmd[i]=='--covdir' | cmd[i]=='-i' ){COV_DIR = cmd[i+1]} # directory where covariates can be found (default: single_cell/project/pop_eQTL/Covariates )
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COV_DIR containing the covariates)
	if (cmd[i]=='--perm' | cmd[i]=='-p' ){perm = cmd[i+1]}
	if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])}
	if (cmd[i]=='--workdir' | cmd[i]=='-w' ){OUT_DIR = cmd[i+1]}
	if (cmd[i]=='--scratch' | cmd[i]=='-h' ){SCRATCH = cmd[i+1]}
	if (cmd[i]=='--run_id' | cmd[i]=='-r' ){RUN_ID = cmd[i+1]}
	if (cmd[i]=='--test' | cmd[i]=='-s' ){TEST = cmd[i+1]}
	if (cmd[i]=='--cellprop' | cmd[i]=='-l' ){ADD_CELLPROP = cmd[i+1]} # cellprop : should cell proportions (computed as mean of the 3 conditions) be added to covariates
	if (cmd[i]=='--expression' | cmd[i]=='-e' ){EXPRESSION_FILE = cmd[i+1]}
	if (cmd[i]=='--metadata' | cmd[i]=='-m' ){META_DATA_FILE = cmd[i+1]}
	if (cmd[i]=='--logfc' | cmd[i]=='-f' ){GET_LOGFC = cmd[i+1]}
}
# myCELLTYPE=gsub(' ','_',myCELLTYPE)

DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/",EVO_IMMUNO_POP_ZEUS)
if(is.null(EXPRESSION_FILE)){
	EXPRESSION_FILE=sprintf('%s/2_population_differences/BatchAdjusted_logCPM_%slibs__per_%s_%s_annotated.tsv.gz',DATA_DIR,NLIBS,CELLTYPE,STATE)
}

if(is.null(COV_DIR)){
	COV_DIR= sprintf("%s/2_population_differences/Covariates",DATA_DIR)
}

###############@ START SCRIPT
options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(rtracklayer))
suppressMessages(library(sva))
suppressMessages(library(ggcorrplot))
suppressMessages(library(tictoc))

####################
# Define RUN NAME
RUN_CELLPROP=ifelse(ADD_CELLPROP,'_19cellPropAdj','')
RUN_LOGFC=ifelse(GET_LOGFC,'_logFC','')
RUN_NAME=sprintf('%s_%s%s_%s%s_%s',CELLTYPE,STATE,RUN_LOGFC,gsub(sprintf('^%s_%s',CELLTYPE,STATE),'',COV_RUN_NAME),RUN_CELLPROP,RUN_ID)

###################
# create output directories
dir.create(sprintf("%s/%s",OUT_DIR,RUN_NAME))
dir.create(sprintf("%s/%s/%s__%s",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE))
dir.create(sprintf("%s/%s/%s__%s/perm",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE))
# dir.create(sprintf('%s/%s/SVs',OUT_DIR,RUN_NAME))

##### define sample meta data
CellMortality=fread(MORTALITY_FILE)
setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
mod=CellMortality[,lm(Mortality~CellCount)]
CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

meta_data=fread(META_DATA_FILE, sep='\t')

##### define IID to use:
MinCell_perCOND=500
keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]

# perm=0
# CHR=22
# myCELLTYPE=5

##### load genotype data
ImputedDIR=sprintf("%s/popCell_data/02_ImputedGenotypes/Imputed/b38/",EVO_IMMUNO_POP_ZEUS)
if(!is.null(CHR)){
	ImputedFILE_CHR=sprintf("Geno_b38_473Ind_3723840snps_chr%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct",CHR)
	VCF=readVcf(sprintf("%s/%s.vcf.gz",ImputedDIR,ImputedFILE_CHR))
	VCF=VCF[!duplicated(rowRanges(VCF)),]
	GenoNum=matrix(unlist(geno(VCF)$DS),nrow(geno(VCF)$DS),ncol(geno(VCF)$DS))
	dimnames(GenoNum)=dimnames(geno(VCF)$GT)
	INFO=as.data.table(info(VCF))
	CHRPOS=rowRanges(VCF)
	# check that there are no multi-allelic SNPs
	ALT_unlist=unlist(rowRanges(VCF)$ALT)
	No_MultiAllelic=( length(ALT_unlist)==nrow(VCF) )
	# if there is no multi-allelic SNP, replace ALT by a DNAStringSet object and unlist AF info field, to speed up & facilitate handling
	if(No_MultiAllelic){
		CHRPOS$ALT=ALT_unlist
		INFO$AF=unlist(INFO$AF)
	}
	CHRPOS=as.data.table(CHRPOS)
	Map=data.table(rsID=rownames(VCF),CHRPOS,INFO)
	Map[,seqnames:=as.character(seqnames)]
	rm(CHRPOS,INFO,ALT_unlist);gc()
	rm(VCF);gc()
}

Map_AF=fread(cmd=sprintf('gunzip -c %s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz | grep -e "ID\\|chr%s"',ImputedDIR,CHR),sep='\t')
Map=merge(Map,Map_AF[,!c('REF','CLASS','ALT','QUAL','FILTER','AF','DR2')],by.x='rsID',by.y='ID')
rm(Map_AF)

##################################################
######### STEP 1:  generate TPM_Matrix  ##########
##################################################

# identify all genotyped individuals from the VCF
indiv_full=colnames(GenoNum)
indiv_full=gsub('(PopCell|EvoImmunoPop)_(ASH|AFB|EUB|PIL|PEL)([0-9]+)(r?)', '\\2\\3', indiv_full)
colnames(GenoNum)=indiv_full

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]
# remove inds with < 500 cells
Expression=Expression[IID%chin%keptIID,]
toc()


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

# select individuals and features to use for eQTL mapping
indiv_to_use=unique(Expression[,.(IID,POP)])
indiv=indiv_to_use$IID

feature_to_use = Expression[,unique(ID)]
toc()

tic('generate feature matrix file and ranktransform it')
#  : Ensembl IDs and Symbol one line per gene, 1 column per individual for the condition/cell_group
Feature_mat=dcast(Expression, ID+Symbol~IID,value.var='logCPM')
setnames(Feature_mat,'ID','feature_id')
set.seed(123)
rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
    }

Feature_mat=t(apply(Feature_mat[match(feature_to_use, Feature_mat$feature_id), mget(indiv)], 1, rankTransform))
rownames(Feature_mat)=feature_to_use
toc()


tic('generating feature_annotation')
#  with EnsemblIDs and Symbol + seqnames, start, end, strand
gtf <- rtracklayer::import(sprintf("%s/single_cell/resources/references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
Feature_annot=as.data.table(gtf)[type=='gene',.(gene_id,gene_name,seqnames, start, end, strand,gene_type)]
Feature_annot[is.na(gene_name),gene_name:=gene_id]
Feature_annot[,seqnames:=as.character(seqnames)]
colnames(Feature_annot)[1]='feature_id'
# for IAV_M and IAV_NS  we have two lines with the same genename (2 transcripts)
# we only keep the longest transcript
Feature_annot=Feature_annot[!(feature_id=="IAV_M" & end==784) & !(feature_id=="IAV_NS" & end==719),]
# we match the remaining features to feature_to_use
Feature_annot=Feature_annot[match(feature_to_use,feature_id),]
toc()

tic('load covariates')
	Covariates=fread(sprintf('%s/%s/Covariates__%s_%s.tsv.gz',COV_DIR,COV_RUN_NAME,myCELLTYPE,mySTATE))
	# Add Cell Mortality
# Covariates=merge(Covariates,CellMortality,by='IID')
#	Covariates=Covariates[,-c('ncells','CellCount')]
	# Add AGGREGATE Cell PROPORTIONS
	# if(ADD_CELLPROP){
	# 	meta=fread(sprintf('%s/single_cell/project/pop_eQTL/data/2_population_differences/sce_clean_metadata.tsv',EVO_IMMUNO_POP_ZEUS))
	# 	CellPct=meta[,.N,by=.(IID,celltype.19level,COND)]
	# 	CellPct=CellPct[,Pct:=N/sum(N),by=.(IID,COND)]
	# 	CellPct_mean=CellPct[,.(Pct=mean(Pct)),by=.(IID,celltype.19level)]
	# 	CellPct_mean=dcast(CellPct_mean,IID~celltype.19level,fill=0)
	# 	Covariates=merge(Covariates,CellPct_mean,by='IID')
	# 	}
		Cov_mat=model.matrix(~0+.,Covariates[,!'IID'])
		rownames(Cov_mat)=Covariates$IID
		Cov_mat=t(Cov_mat[indiv,])
toc()

########################################################################
##										Genes by Pop - Standard - Cis 									##
########################################################################

RESCIS=list()
QTL_pvalues=list()
tim=Sys.time()

wSNP=which(Map$max_MAF>minFreq)#[1:10000]
if(TEST){wSNP=wSNP[1:10000]}

CPM_Mat = SlicedData$new()
CPM_Mat$CreateFromMatrix(Feature_mat)
CPM_Mat$ResliceCombined(sliceSize = 5000)
show(CPM_Mat)

Cvrt=SlicedData$new()
Cvrt$CreateFromMatrix(Cov_mat)

if(perm>0){
		set.seed(perm)
		indiv[grep('AFB',indiv)]=sample(indiv[grep('AFB',indiv)])
		indiv[grep('EUB',indiv)]=sample(indiv[grep('EUB',indiv)])
		indiv[grep('ASH',indiv)]=sample(indiv[grep('ASH',indiv)])
		}

GenoNum_Mat = SlicedData$new()
GenoNum_Mat$CreateFromMatrix(GenoNum[wSNP,match(indiv,colnames(GenoNum))])
GenoNum_Mat$ResliceCombined(sliceSize = 5000)

res=Matrix_eQTL_main(GenoNum_Mat,
	                     CPM_Mat,
	                     Cvrt,
	                     output_file_name=NULL,
	                     pvOutputThreshold=0,
	                     output_file_name.cis=NULL,
	                     pvOutputThreshold.cis=pvCis,
	                     cisDist=CIS_DIST,
	                     snpspos=Map[wSNP,c("rsID","seqnames","start")],
	                     genepos=Feature_annot[,c('feature_id','seqnames','start','end')],
	                     min.pv.by.genesnp=TRUE)

	if(nrow(res$cis$eqtls)>0){
		RESCIS=res$cis$eqtls
		RESCIS$celltype=myCELLTYPE
	}
    Pvalues=merge(data.table(feature_id=names(res$cis$min.pv.gene),
                                    pvalue=res$cis$min.pv.gene,
                                    celltype=myCELLTYPE),Feature_annot)
    QTL_pvalues=Pvalues[seqnames==paste('chr',CHR,sep=''),mget(c('feature_id','pvalue','celltype'))]
	print(Sys.time()-tim)

#QTL_pvalues=dcast(QTL_pvalues,feature_id+QTL_type~condition,value.var='pvalue')
#QTL_pvalues[,Min_Pvalue:=min(pvalue_NS,pvalue_LPS,pvalue_PAM3CSK4,pvalue_R848,pvalue_IAV))
fwrite(QTL_pvalues,file=sprintf("%s/%s/%s__%s/perm/eQTL_ALL_chr%s_bestP_perFeature_perm%s.txt.gz",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE,CHR,perm),sep='\t')

if(perm<=1){
    RESCIS$snps=as.character(RESCIS$snps)
    RESCIS$gene=as.character(RESCIS$gene)
		dfFull=res$param$dfFull
		RESCIS$R2=RESCIS$statistic^2/(RESCIS$statistic^2+dfFull)
		RESCIS$se=RESCIS$beta/RESCIS$statistic
		RESCIS$FDR=NULL
    wSNP=match(RESCIS$snp,Map$rsID)
    wFeature=match(RESCIS$gene,Feature_annot$feature_id)
    RESCIS$CisDist=ifelse(Map[wSNP,seqnames]!=Feature_annot[wFeature,seqnames],Inf,
    		ifelse(Map[wSNP,end] < Feature_annot[wFeature,start],
    															Feature_annot[wFeature,start]-Map[wSNP,end],
											 						ifelse(Feature_annot[wFeature,end]<Map[wSNP,start],
																	Map[wSNP,start]-Feature_annot[wFeature,end],0))
    		)
	if(perm==0){
			fwrite(RESCIS,file=sprintf("%s/%s/%s__%s/eQTL_ALL_chr%s_assoc.txt.gz",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE,CHR),sep='\t')
		}else{
			fwrite(RESCIS,file=sprintf("%s/%s/%s__%s/eQTL_ALL_chr%s_assoc_permuted.txt.gz",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE,CHR),sep='\t')
		}
	}

q('no')
