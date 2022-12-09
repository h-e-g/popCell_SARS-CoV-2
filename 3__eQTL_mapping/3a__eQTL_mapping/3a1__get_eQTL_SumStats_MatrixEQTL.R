################################################################################
################################################################################
# File name: 3a1__get_eQTL_SumStats_MatrixEQTL.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use MatrixEQTL models to estimate eQTL effect sizes, R2, and p-value
# same script is used for reQTLs
# Effector script
################################################################################
################################################################################

################################################################################
# Setup
##################### Expected inputs:
# Expression file: with logCPM of all genes, per IID x celltype x state
# metadata file: with 1 line per IID x COND x RUN, used to annotate samples

# define default values
TEST=FALSE # if true, run on a samll susbet of SNPs to test code
CELLTYPE='lineage'
STATE='condition'
COV_DIR=NULL

DATA_DIR = "3__eQTL_mapping"
OUT_DIR = "3_eQTL_mapping/sumStats"
COV_RUN_NAME='lineage_condition__CellPropLineage_SVs'
GET_LOGFC=FALSE
EXPRESSION_FILE=NULL
RUN_ID="220409" # <<<<update when rerunning to generate files in a separate folder >>>>


CIS_DIST=1e6 # Distance in Cis to consider in the analysis.
pvCis=1 # min P-value in cis to report a variant
minFreq=0.05 # min allele frequency (per population) to consider a variant
NLIBS=125 # number of libraries considered to the run
TEST=FALSE #for testing purpose. Set to true to run on the first 100 SNPs


cmd=commandArgs()
for (i in 1:length(cmd)){
	if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use to define pseudoBulk. Will be used for naming of output files
	if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]}
	if (cmd[i]=='--mystate' | cmd[i]=='-o' ){mySTATE = cmd[i+1]} # condition for which the eQTL mapping will be done
	if (cmd[i]=='--mycelltype' | cmd[i]=='-y' ){myCELLTYPE = cmd[i+1]} # celltype for which the eQTL mapping will be done (should be one level of CELLTYPE variable)
	if (cmd[i]=='--covname' | cmd[i]=='-v' ){COV_RUN_NAME = cmd[i+1]} # name of the set of covariates to be used (COV_RUN_NAME) (subdirectory of COVAR_DIR/Covariates containing the desired covariates)
	if (cmd[i]=='--perm' | cmd[i]=='-p' ){perm = cmd[i+1]} # should the data be permuted ? 0 for observed data. if >0, the value will be used as seed for the permutation
	if (cmd[i]=='--run_id' | cmd[i]=='-r' ){RUN_ID = cmd[i+1]}
	if (cmd[i]=='--logfc' | cmd[i]=='-f' ){GET_LOGFC = cmd[i+1]} # should the mapping be done on the response ?
}

if(is.null(EXPRESSION_FILE)){
	EXPRESSION_FILE=sprintf('%s/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz',EXPR_DIR,NLIBS,CELLTYPE,STATE)
}

COV_DIR= sprintf("%s/Covariates/",COVAR_DIR)


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
RUN_LOGFC=ifelse(GET_LOGFC,'_logFC','')
RUN_NAME=sprintf('%s_%s%s_%s_%s',CELLTYPE,STATE,RUN_LOGFC,gsub(sprintf('^%s_%s',CELLTYPE,STATE),'',COV_RUN_NAME),RUN_ID)

###################
# create output directories
dir.create(OUT_DIR)
dir.create(sprintf("%s/%s",OUT_DIR,RUN_NAME))
dir.create(sprintf("%s/%s/%s__%s",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE))
dir.create(sprintf("%s/%s/%s__%s/perm",OUT_DIR,RUN_NAME,myCELLTYPE,mySTATE))

##### load genotype data
if(!is.null(CHR)){
	ImputedFILE_CHR=sprintf("Geno_b38_473Ind_3723840snps_chr%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct",CHR)
	VCF=readVcf(sprintf("%s/%s.vcf.gz",Genotype_DIR,ImputedFILE_CHR))
	VCF=VCF[!duplicated(rowRanges(VCF)),]
	GenoNum=matrix(unlist(geno(VCF)$DS),nrow(geno(VCF)$DS),ncol(geno(VCF)$DS))
	dimnames(GenoNum)=dimnames(geno(VCF)$GT)
	# create a genotype annotation map
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
# read file containing MAF for each population
Map_AF=fread(cmd=sprintf('%s/Map_allCHR_b38_AF_3pops.tsv.gz',Genotype_DIR,CHR),sep='\t')
Map=merge(Map,Map_AF[,.(ID, maf_AFB,maf_EUB,maf_ASH, max_MAF=pmax(maf_AFB,maf_EUB,maf_ASH))],by.x='rsID',by.y='ID')
rm(Map_AF)

##################################################
######### STEP 1:  generate TPM_Matrix  ##########
##################################################

# identify all genotyped individuals from the VCF
indiv_full=colnames(GenoNum)
indiv_full=gsub('(PopCell|EvoImmunoPop)_(ASH|AFB|EUB)([0-9]+)(r?)', '\\2\\3', indiv_full)
colnames(GenoNum)=indiv_full

tic('loading batch-adjusted CPM data')
Expression=fread(file=EXPRESSION_FILE)
Expression=Expression[,-c('ncells','Age','Gender')]
# remove inds with < 500 cells
remove_IID=fread('1__transcriptome_processing/data/low_cellcount_donors.tsv',header=FALSE)$V1
Expression=Expression[!IID%chin%remove_IID,]
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
gtf <- rtracklayer::import("references/human_iav_sars/genes/genes.gtf")
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
