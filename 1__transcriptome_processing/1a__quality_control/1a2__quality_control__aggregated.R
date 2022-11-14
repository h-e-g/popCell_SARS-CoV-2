################################################################################
################################################################################
# File name: 1a2__quality_control__aggregated.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Library-aggregated transcriptome quality control
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("./1a__quality_control__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="MISC"
source(sprintf("%s/shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_functions.R",MISC_DIR))

# read-in library ID
args <- commandArgs(TRUE)
LIB=args[1]

################################################################################
# Data loading

print("Loading filtered data for all libraries.")

# check which libraries have been processed and how many cells they contain
processedLIBs=dir(sprintf('%s/',LIB_QC_DIR),pattern='cells_sce.RDS')

# remove libraries from time points 0 and 24
processedLIBs=processedLIBs[!grepl("L113|L114|L121|L122|L123|L124|L125|L126",processedLIBs)]

# remove libraries from with mixed-condition samples
processedLIBs=processedLIBs[!grepl("L117|LNA",processedLIBs)]

LIB=gsub('(L[0-9]+)_([0-9]+)cells_sce.RDS','\\1',processedLIBs)
nCELL=as.numeric(gsub('(L[0-9]+)_([0-9]+)cells_sce.RDS','\\2',processedLIBs))
nLIBs=length(LIB)

# load library metadata
meta_data=fread("../../0__barcode_processing/data/experimental_metadata.tsv")

# load all SingleCellExperiment objects
sce <- lapply(1:length(processedLIBs),function(l){
  readRDS(sprintf('%s/%s_%scells_sce.RDS',LIB_QC_DIR,LIB[l],nCELL[l]))
})

################################################################################
# Data aggregation

# rescale library size factors to account for sequencing depth differences across batches
rescaled <- multiBatchNorm(sce)

# model gene variance per batch to account for different mean-variance relations and combine results
all.pois <- lapply(rescaled, modelGeneVarByPoisson)
combined.pois<- combineVar(all.pois)

# extract highliy variables genes
chosen.hvgs <- combined.pois$bio > 0.001 | combined.pois$mean>0.01

for( i in 1:length(rescaled) ){
  rescaled[[i]]$batch=names(rescaled)[i]
}

# aggregate rescaled SingleCellExperiment objects and log-normalize rescaled counts
rescaled=CombineSCE(rescaled)
rescaled <- logNormCounts(rescaled)

# compute denoised PCs from highly variable genes
rescaled=denoisePCA(rescaled, technical=combined.pois, subset.row=chosen.hvgs,max.rank=100)

# filter out non-highly variable genes
rescaled=rescaled[chosen.hvgs,]

# compute UMAP on rescaled and aggregated data
rescaled <- runUMAP(rescaled, dimred="PCA")

# set batch labels

meta_data <- as.data.table(rescaled@colData)
meta_data[,run:=case_when(
  LIB%chin%paste0("L",1:8)~"1",
  LIB%chin%paste0("L",9:16)~"2",
  LIB%chin%paste0("L",17:24)~"3",
  LIB%chin%paste0("L",25:32)~"4",
  LIB%chin%paste0("L",33:40)~"5",
  LIB%chin%paste0("L",41:48)~"6",
  LIB%chin%paste0("L",49:56)~"7",
  LIB%chin%paste0("L",57:64)~"8",
  LIB%chin%paste0("L",65:72)~"9",
  LIB%chin%paste0("L",73:80)~"10",
  LIB%chin%paste0("L",81:88)~"11",
  LIB%chin%paste0("L",89:96)~"12",
  LIB%chin%paste0("L",97:104)~"13",
  LIB%chin%paste0("L",105:112)~"14",
  LIB%chin%paste0("L",113:126)~"TC",
  LIB%chin%paste0("L",127:130)~"R1",
  LIB%chin%paste0("L",131:134)~"R2"
)]
rescaled$run <- meta_data[,setNames(run,Barcode)][rescaled$Barcode]

# use Harmony to remove variation associated to experimental run from denoised PCs
reducedDim(rescaled, "harmony_PCA_run") <-  HarmonyMatrix(reducedDim(rescaled,"PCA"),rescaled$run,"run",do_pca=FALSE)

# find clusters on Harmony-corrected PC space
clusters <- custom_clustering(reducedDim(rescaled,"harmony_PCA_run"),25)
rescaled$cluster_seurat_harmony=clusters

# find clusters on uncorrected PC space
clusters <- custom_clustering(reducedDim(rescaled,"PCA"),25)
rescaled$cluster_seurat=clusters

# compute UMAP from Harmony-corrected PC space
rescaled <- runUMAP(rescaled, dimred="harmony_PCA_run",name="UMAP_harmony_run")

################################################################################
# Write results

AGGR_QC_FILE="sce_rescaled_full_125libs_updated__harmony_run.RDS"

saveRDS(rescaled,sprintf("%s/data/%s",AGGR_QC_DIR,AGGR_QC_FILE))
