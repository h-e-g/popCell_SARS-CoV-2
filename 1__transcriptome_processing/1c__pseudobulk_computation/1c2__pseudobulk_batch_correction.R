################################################################################
################################################################################
# File name: 1c2__pseudobulk_batch_correction.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Pseudobulk batch correction
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./1c__pseudobulk_computation.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))

# read-in library ID
args <- commandArgs(TRUE)
LIB=args[1]

# define default parameters
NLIBS=125
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition'
ncell_threshold=1
OUT_DIR = "1__transcriptome_processing" # output directory

# update parameter values based on arguments provided in 1c2__pseudobulk_batch_correction.sh
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
        if (cmd[i]=='--meta' | cmd[i]=='-m' ){META_DATA_FILE = cmd[i+1]} # path to meta data file with one line per IID, COND, LIB
        if (cmd[i]=='--count' | cmd[i]=='-c' ){COUNT_DATA_FILE = cmd[i+1]} # path to count & CPM data file withs one line per LIB, IID, CELLTYPE, COND, and gene.
  if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use. Will be used for naming of output files (either lineage or celltype)
  if (cmd[i]=='--nlibs' | cmd[i]=='-n' ){NLIBS = cmd[i+1]} # number of libraries in the dataset (for naming)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output directory
}

################################################################################
# Data loading

# load metadata file: one line per IID, COND, LIB
META_DATA_FILE="1__transcriptome_processing/data/experimental_metadata.tsv"

# load pseudobulk data created in 1c1__pseudobulk_computation.R : one line per LIB, IID, CELLTYPE, STATE, and gene
COUNT_DATA_FILE=sprintf('1__transcriptome_processing/data/pseudobulk_%slibs__per_%s_%s_IID_and_LIB.tsv.gz',NLIBS,CELLTYPE,STATE)

# load library meta_data
meta_data=fread(META_DATA_FILE,sep='\t')

# load pseudobulk count data
Counts=fread(COUNT_DATA_FILE)
Counts=merge(Counts,meta_data[LIB%in%unique(Counts$LIB),.(LIB,IID,COND,POP,p3_date,replicate,flowcellID)],by=c('LIB','IID'),allow.cartesian=T)
Counts[,p3_date:=as.factor(p3_date)]

# keep reasonably expressed genes (meanCPM>1 in at least one condtion/celltype) to avoid features with too few reads
# compute mean CPM per (gene by cell type by condition)
MeanCPM=Counts[,.(meanCPM=mean(CPM)),by=.(ID,Symbol,celltype, state)]
expressed_genes=MeanCPM[meanCPM>1,.N,by=ID][,ID]

################################################################################
# Mixed linear model to estimate and remove batch effects

# declare MLM function
estimate_donor_effect=function(CPM, IID, LIB, p3_date, flowcellID, logged=FALSE, addNoise=0, test_ID=''){
  # CPM: count per million
  # IID: indiividual ID
  # LIB: library ID
  # p3_date: experiment batch (RUN)
  # flowcellID: flowcellID (eg. HGNTGCCX2)
  # logged: should CPM be log transformed (with an offset of 1)
  # addNoise : if>0, a small amount of noise (with variance addNoise) will be added to reduce the risk of numerical issues
  # testID : character string that will be printed for each call (allowing to monitor progress)

  cat(unique(test_ID),'\n')
  flush.console()

  if(logged){
    transformedCPM=log2(1+CPM)
  }else{
    transformedCPM=CPM
  }

  if(addNoise>0){
    transformedCPM=transformedCPM+rnorm(length(transformedCPM),0,addNoise)
  }

  model=try(lmer(transformedCPM~(1|IID)+(1|LIB)+(1|p3_date)+(1|flowcellID)))
  if(class(model)=='try-error'){
      # if something goes wrong (e.g., 0s in one cell type x cell state)
      # set estimated R2s to negative values
      R2_vars=rep(-999,5)
      names(R2_vars)=c('IID','LIB','p3_date','flowcellID','residual')
      # return untrandformed data
      uIID=unique(IID)
      DT=data.table(residuals=transformedCPM,IID)
      DT=DT[,.(residuals=mean(residuals)),by=IID]
      IID_effect=DT[match(uIID,DT$IID),residuals]
      names(IID_effect)=uIID
  }else{
    # compute model summary & define R2
    Model_summary=try(summary(model))
    if(class(Model_summary)=='try-error'){
        R2_vars=rep(-9999,5)
        names(R2_vars)=c('IID','LIB','p3_date','flowcellID','residual')
      }else{
        se_vars=c(sqrt(unlist(Model_summary$varcor)),residual=summary(model)$sigma)
        R2_vars=t(se_vars^2)/apply(t(se_vars)^2,1,sum)
        names(R2_vars)=c('IID','LIB','p3_date','flowcellID','residual')
    }
    # extract model residuals
    DT=data.table(residuals=residuals(model),IID)
    # average residuals by IID
    DT=DT[,.(residuals=mean(residuals)),by=IID]
    uIID=unique(IID)
    # extract intercept + IID effects +averaged residuals
    IID_effect=model@beta[1] + ranef(model)$IID[uIID,1] + DT[match(uIID,DT$IID),residuals]
    names(IID_effect)=uIID
    # add residuals
  }
  # return R2 and IID effects as a list (variable & value)
  list(value=c(R2_vars,IID_effect),variable=c(names(R2_vars),names(IID_effect)))
}

# apply MLM function
Prov=Counts[ID%in%expressed_genes,estimate_donor_effect(CPM, IID, LIB, p3_date, flowcellID, logged=TRUE, test_ID=paste(celltype, state ,ID ,Symbol,sep='_')),by=.(celltype,state,ID,Symbol)]

# extract variance estimates and save them to disk
Var_estimates=Prov[variable%chin%c('IID','LIB','p3_date','flowcellID','residual'),]
fwrite(Var_estimates,file=sprintf('%s/data/varexplained_by_IID_and_batch_%slibs__per_%s_%s.tsv.gz',OUT_DIR,NLIBS,CELLTYPE,STATE),sep='\t')

# extract donor effects and save them to disk
BatchAdjusted_logCPM=Prov[!variable%chin%c('IID','LIB','p3_date','flowcellID','residual'),]
setnames(BatchAdjusted_logCPM,'variable','IID')
setnames(BatchAdjusted_logCPM,'value','logCPM')
# setnames(IID_effect,CELLTYPE_lowercase,'celltype')

tic('adding cell count information and meta_data to main table')
#### compute cellcounts per cluster/IID
Cells_per_sample=Counts[,.(ncells=sum(ncells)),by=.(celltype,state,IID,ID)]

#### add cellcounts to main table
BatchAdjusted_logCPM=merge(BatchAdjusted_logCPM, Cells_per_sample,by=c('celltype','state','ID','IID'))
#### add meta_data to main table
BatchAdjusted_logCPM=merge(BatchAdjusted_logCPM, meta_data[!duplicated(IID),.(IID,Age,Gender,POP)],by=c('IID'))

fwrite(BatchAdjusted_logCPM,file=sprintf('%s/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz',EXPR_DIR,NLIBS,CELLTYPE,STATE),sep='\t')
toc()
