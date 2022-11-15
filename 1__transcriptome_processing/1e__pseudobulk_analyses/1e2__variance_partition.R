################################################################################
################################################################################
# File name: 1e2__variance_partition.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute CAR scores for cell type and condition
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./1a__quality_control__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))

# define default parameters
NLIBS=125
CELLTYPE='lineage' # celltype variable to use. Will be used for naming of output files
STATE='condition'

# load batch-adjusted counts computed in 1c2__pseudobulk_batch_correction.R
Expr=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))

################################################################################
# Fraction of gene expression variation explained by condition and cell type

# build covariate matrices
celltype_covariates=dcast(
  Expr[Symbol=="FGR",.(IID,celltype,state,dummy=1)],
  IID+celltype+state~celltype,
  value.var="dummy",fill=0
)[,-"T.CD4"]

condition_covariates=dcast(
  Expr[Symbol=="FGR",.(IID,celltype,state,dummy=1)],
  IID+celltype+state~celltype+state,
  value.var="dummy",fill=0
)
condition_covariates=condition_covariates[,mget(colnames(condition_covariates)[!grepl("NS",colnames(condition_covariates))])]

response_vectors_y=dcast(
  Expr[,.(IID,celltype,state,Symbol,logCPM)],
  IID+celltype+state~Symbol,
  value.var="logCPM"
)

pseudobulk_model <- list(
      cbind(as.matrix(celltype_covariates[,-c("IID","celltype","state")]),
            as.matrix(condition_covariates[,-c("IID","celltype","state")])),
      as.matrix(response_vectors_y[,-c("IID","celltype","state")])
)
names(pseudobulk_model) <- c("covariates_matrix_x","response_vectors_y")

# compute CAR scores
car <- lapply(X=1:ncol(pseudobulk_model[["response_vectors_y"]]),FUN=function(g){
    print(g)
    carscore(
      pseudobulk_model[["covariates_matrix_x"]],
      pseudobulk_model[["response_vectors_y"]][,g],
      lambda=0,
      verbose=F
    )
})
names(car) <- colnames(pseudobulk_model[["response_vectors_y"]])

# reshape and format data
var_explained_celltype <-
  sapply(X=1:length(car),FUN=function(g){sum(car[[g]][1:4]**2)})
names(var_explained_celltype) <- names(car)
var_explained_condition <-
  sapply(X=1:length(car),FUN=function(g){sum(car[[g]][5:14]**2)})
names(var_explained_condition) <- names(car)

variance_explained=data.table(Symbol=names(var_explained_celltype),celltype=var_explained_celltype,condition=var_explained_condition)
variance_explained=melt(variance_explained,measure.vars=2:3,variable.name="covariate",value.name="variance_explained")

fwrite(variance_explained,"1__transcriptome_processing/data/variance_explained_celltype_condition.tsv.gz",sep="\t")
