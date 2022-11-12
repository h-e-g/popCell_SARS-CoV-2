######## objective of the script:
# define expressed genes (mean CPM>1 in at least one CELLTYPE & STATE)
# for each level of CELLTYPE and STATE, compute logTPM and correct them for experimental BATCH / library effects
# add meta data (age, Gender, cell number) to the resulting table
# estimate variance explained by batch/library/IID...


###### requirements:
### define CELLTYPE and STATE
# ideally CELLTYPE and CELLSTATE, should not contain any '_'
# eg. CELLTYPE='celltype.19level'
#     STATE='activation.state'
#
# or CELLTYPE='celltype.8level'
#     STATE='condition'

### define path to the count object generated in script 6a

### define path to an up-to-date metadata file

### set NLIBS to the number of libraries in the current iteration (for naming purposes)

### set OUT_DIR to the desired outpout directory ( evo_immuno_pop/single_cell/project/pop_eQTL/data by default)

.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggrastr))
suppressMessages(library(scales))
suppressMessages(library(data.table))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DropletUtils))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(BiocNeighbors))
suppressMessages(library(SoupX))
suppressMessages(library(CelliD))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(batchelor))
suppressMessages(library(kBET))
suppressMessages(library(lme4))
suppressMessages(library(tictoc))

EVO_IMMUNO_POP_ZEUS = "/pasteur/zeus/projets/p02/evo_immuno_pop"
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences",EVO_IMMUNO_POP_ZEUS)
OUT_DIR = DATA_DIR # output directory (default: DATA_DIR)
NLIBS=125 # number of libraries in the dataset (for naming)
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition' # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
ncell_threshold=1
# internally these variables are stored in celltype and state

# update parameter values based on provided arguments
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--meta' | cmd[i]=='-m' ){META_DATA_FILE = cmd[i+1]} # path to meta data file with one line per IID, COND, LIB
	if (cmd[i]=='--count' | cmd[i]=='-c' ){COUNT_DATA_FILE = cmd[i+1]} # path to count & CPM data file withs one line per LIB, IID, CELLTYPE, STATE, and gene.
  if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # celltype variable to use. Will be used for naming of output files
	if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
  if (cmd[i]=='--nlibs' | cmd[i]=='-n' ){NLIBS = cmd[i+1]} # number of libraries in the dataset (for naming)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output directory
}

# metadata file: one line per IID, COND, LIB
META_DATA_FILE=sprintf('%s/popCell_data/00_CRF/scrnaseq_LIB_IID_COND_based_metadata_full_long_withQC_v2.tsv',EVO_IMMUNO_POP_ZEUS)

# count & CPM data : one line per LIB, IID, CELLTYPE, STATE, and gene.
COUNT_DATA_FILE=sprintf('%s/Counts_and_CPM__pseudoBulk_full_%slibs__per_%s_%s_IID_and_LIB.tsv.gz',DATA_DIR,NLIBS,CELLTYPE,STATE)

# load library meta_data generated in previous script (4a)
meta_data=fread(META_DATA_FILE,sep='\t')

# load count data generated in previous script (6a)

Counts=fread(COUNT_DATA_FILE)
Counts=merge(Counts,meta_data[LIB%in%unique(Counts$LIB),.(LIB,IID,COND,POP,p3_date,replicate,flowcellID)],by=c('LIB','IID'),allow.cartesian=T)
Counts[,p3_date:=as.factor(p3_date)]

######## Identify reasonnably expressed genes to avoid genes with too few reads...
#### (meanCPM>1 in at least one condtion/celltype)

# compute mean CPM per (gene x  cell type x cell state)
MeanCPM=Counts[,.(meanCPM=mean(CPM)),by=.(ID,Symbol,celltype, state)]
expressed_genes=MeanCPM[meanCPM>1,.N,by=ID][,ID]

# write down this list
fwrite(data.table(expressed_genes),file=sprintf('%s/expressed_genes_CPMover1_%slibs__per_%s_and_%s.tsv',OUT_DIR, NLIBS, CELLTYPE, STATE),sep='\t')
# TODO: create another file excluding chr X and Y genes to define the list of genes to use in next steps


##########################################################################
##########################################################################
########  use mixed model to estimate and remove batch effects  ##########
##########################################################################
##########################################################################

explainVar_estimate_IIDeffect=function(CPM, IID, LIB, p3_date, flowcellID, logged=FALSE, addNoise=0, test_ID=''){
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
      ## if something goes wrong (eg. 0s in one cell type x cell state)
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


# apply explainVar_estimate_IIDeffect
tim=Sys.time()
Prov=Counts[ID%in%expressed_genes,explainVar_estimate_IIDeffect(CPM, IID, LIB, p3_date, flowcellID, logged=TRUE, test_ID=paste(celltype, state ,ID ,Symbol,sep='_')),by=.(celltype,state,ID,Symbol)]
print(Sys.time()-tim) # 18 hours for 12845 genes


# extract Variance estimates and save them to disk
Var_estimates=Prov[variable%chin%c('IID','LIB','p3_date','flowcellID','residual'),]
fwrite(Var_estimates,file=sprintf('%s/VarianceExplained_by_IID_and_batch_%slibs__per_%s_%s.tsv.gz',OUT_DIR,NLIBS,CELLTYPE,STATE),sep='\t')

# extract IID effects and save them to disk
BatchAdjusted_logCPM=Prov[!variable%chin%c('IID','LIB','p3_date','flowcellID','residual'),]
setnames(BatchAdjusted_logCPM,'variable','IID')
setnames(BatchAdjusted_logCPM,'value','logCPM')
# setnames(IID_effect,CELLTYPE_lowercase,'celltype')

fwrite(BatchAdjusted_logCPM,file=sprintf('%s/BatchAdjusted_logCPM_%slibs__per_%s_%s.tsv.gz',OUT_DIR,NLIBS,CELLTYPE,STATE),sep='\t')


tic('adding cell count information and meta_data to main table')
#### compute cellcounts per cluster/IID
Cells_per_sample=Counts[,.(ncells=sum(ncells)),by=.(celltype,state,IID,ID)]

#### add cellcounts to main table
BatchAdjusted_logCPM=merge(BatchAdjusted_logCPM, Cells_per_sample,by=c('celltype','state','ID','IID'))
#### add meta_data to main table
BatchAdjusted_logCPM=merge(BatchAdjusted_logCPM, meta_data[!duplicated(IID),.(IID,Age,Gender,POP)],by=c('IID'))

fwrite(BatchAdjusted_logCPM,file=sprintf('%s/BatchAdjusted_logCPM_%slibs__per_%s_%s_annotated.tsv.gz',OUT_DIR,NLIBS,CELLTYPE,STATE),sep='\t')
toc()
