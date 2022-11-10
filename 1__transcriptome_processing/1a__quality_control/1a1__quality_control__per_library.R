################################################################################
################################################################################
# File name: 1a1__quality_control__per_library.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-library transcriptome quality control 
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

# read-in library ID 
args <- commandArgs(TRUE)
LIB=args[1]

################################################################################
# Per-library quality control: round 1

print(sprintf("Running first QC on L%s.", LIB))

# read-in the data and compute filtering metrics

print("Reading STARsolo outputs.")

sce <- sprintf("%s/%s.Solo.out/Gene/filtered/",ALIGN_DIR,LIB)
sce <- read10xCounts(sce,col.names=TRUE)
print(sce)

sce$LIB <- LIB

print("Identifying gene subsets to be used for QC metrics.")

is.MT <- grepl("^MT-", rowData(sce)$Symbol)
print("Number of mitochondrial transcripts:")
table(is.MT)

is.IAV COVrepl("^IAV_", rowData(sce)$Symbol)
print("Number of IAV transcripts:")
table(is.IAV)

is.COV <- grepl("^SARS_", rowData(sce)$Symbol)
print("Number of COV transcripts:")
table(is.COV)

# calculate QC metrics

print("Computing QC metrics.")

sce <- addPerCellQC(sce,subsets=list(MT=is.MT,IAV=is.IAV,COV=is.COV))

# flag library level QC (soft filters) and apply hard filters (library-independent)

# soft filters (SF)

countsSF <- quantile(sce$total, probs=c(0.01, 0.05, 0.95, 0.99))
print("Quantiles for total counts:")
print(countsSF)

genesSF <- quantile(sce$detected, probs=c(0.01, 0.05, 0.95, 0.99))
print("Quantiles for gene features identified:")
print(genesSF)

mitoSF <- quantile(sce$subsets_MT_percent, probs=c(0.01, 0.05, 0.95, 0.99))
print("Quantiles for mitochondrial content:")
print(mitoSF)

# flag library-level 1% and 5% QC outliers
sce$lib_QC1 <- ifelse(sce$total >= countsSF["1%"] & sce$detected >= genesSF["1%"] & sce$subsets_MT_percent >= mitoSF["99%"], "Pass", "Fail")
sce$lib_QC5 <- ifelse(sce$total >= countsSF["5%"] & sce$detected >= genesSF["5%"] & sce$subsets_MT_percent >= mitoSF["95%"], "Pass", "Fail")

# hard filters (HF)
countsHF=1500
genesHF=500
mitoHF=20

sce$QC <- ifelse(sce$total >= countsHF & sce$detected >= genesHF & sce$subsets_MT_percent <= mitoHF, "Pass", "Fail")
print(sprintf("Filtering cells with less than %s total counts, less than %s genes, and/or higher that %s percent mitochondrial content"),countsHF,genesHF,mitoHF)

print("Filtered cell counts:")
print(table(sce$QC))

sce <- sce[,sce$QC=="PASS"]

print("Filtered data set:")
print(sce)

################################################################################
# Post-filtering normalization: round 1

# compute size factors and log-normalize

print("Computing size factors.")

set.seed(1000)
quickclusters <- quickCluster(sce,min.mean=0.1,method="igraph")

sce$quickclusters <- as.factor(paste(gsub("^","L",LIB),quickclusters,sep="_"))
sce <- computeSumFactors(sce, cluster=quickclusters,min.mean=0.1)

print("Summary of deconvoluted size factors:")
summary(sizeFactors(sce))

print("Log-normalizing.")

sce <- logNormCounts(sce)

################################################################################
# Feature selection: round 1

print("Decomposing gene log-expression variance into technical and biological components.")

set.seed(1000)
dec.pois.sce <- modelGeneVarByPoisson(sce)

print("Extracting highly-variable genes.")

hvg.sce.var <- getTopHVGs(dec.pois.sce,fdr.threshold=0.01)

print("Genes with biological variance (FDR 1%):")
print(length(hvg.sce.var))

################################################################################
# Dimensionality reduction: round 1

# principal components analysis and log-expression data denoising by removing first PCs

print("Performing library-level PCA.")

set.seed(1000)
sce <- denoisePCA(sce, technical=dec.pois.sce, subset.row=hvg.sce.var)

################################################################################
# Cryptic doublet removal

# read-in Demuxlet results and add calls to SingleCellExperiment object

print("Adding Demuxlet's calls.")

DMX <- data.frame(fread(sprintf('%s/demuxlet/%s/%s_demuxlet_fullGenotypes.best',DMX_DIR,LIB,LIB)))
DMX$SNG.BEST.GUESS=sapply(strsplit(as.character(DMx$SNG.BEST.GUESS),"_"),"[[",2)
DMX$IID <- ifelse(
  DMX$DROPLET.TYPE=="SNG",DMX$SNG.BEST.GUESS,ifelse(
    DMX$DROPLET.TYPE=="DBL","doublet","ambiguous"
  )
)

sce$DROPLET.TYPE <- DMX[match(sce$Barcode,DMx$BARCODE),'DROPLET.TYPE']
sce$IID <- DMX[match(sce$Barcode, DMx$BARCODE), 'IID']

sce <- sce[,which(!is.na(sce$DROPLET.TYPE))]

# identify and remove cryptic doublets

print("Performing nearest-neighbor search for doublets arising from the same donor.")

set.seed(1000)

nn <- findKNN(reducedDim(sce,"PCA"),k=25,BNPARAM=KmknnParam())

sce$DROPLET.TYPE <- factor(sce$DROPLET.TYPE,levels=c("SNG", "DBL", "AMB"))
meta <- data.frame(colData(sce))

# count number of doublets in the 25-nearest neighbors of each barcode in PCA space

dbl=c()
for (i in 1:dim(sce)[2]) {
  dbl=c(dbl, unname(table(meta[nn$index[i,],'DROPLET.TYPE'])[2]))
}

sce$nn_dbl <- dbl

print("Removing non-singlet and suspected cryptic doublet barcodes.")

clean <- sce[,sce$DROPLET.TYPE=="SNG" & sce$nn_dbl<=5]

################################################################################
# Post-filtering normalization: round 2

# compute size factors and log-normalize

print("Computing size factors on clean data.")

set.seed(1000)
quickclusters <- quickCluster(clean,min.mean=0.1,method="igraph")

clean$quickclusters <- as.factor(paste(gsub("^","L",LIB),quickclusters,sep="_"))
clean <- computeSumFactors(clean,cluster=quickclusters,min.mean=0.1)

print("Summary of deconvoluted size factors on clean data:")
summary(sizeFactors(clean))

print("Log-normalizing clean data.")

clean <- logNormCounts(clean)

################################################################################
# Feature selection: round 2

print("Decomposing gene log-expression variance into technical and biological components.")

set.seed(1000)
dec.pois.clean <- modelGeneVarByPoisson(clean)

print("Extracting highly-variable genes.")

hvg.clean.var <- getTopHVGs(dec.pois.clean,fdr.threshold=0.01)

print("Genes with biological variance (FDR 1%):")
print(length(hvg.clean.var))

################################################################################
# Dimensionality reduction: round 2

# principal components analysis and log-expression data denoising by removing first PCs

print("Performing library-level PCA.")

set.seed(1000)
clean <- denoisePCA(clean, technical=dec.pois.clean, subset.row=hvg.clean.var)

# use Harmony to reduce variation related to individual and condition, compute t-SNE and UMAP

print("Harmony adjustment of denoised PCs by donor.")

set.seed(1000)
reducedDim(clean,"harmony_PCA") <- HarmonyMatrix(
  reducedDim(clean,"PCA"),clean$IID,"IID",
  do_pca=FALSE
)

print(str(clean))

print("Computing library-level t-SNE on cleaned data.")

set.seed(1000)
clean <- runTSNE(clean,dimred="PCA")

print("Computing library-level UMAP on cleaned data.")

set.seed(1000)
clean <- runUMAP(clean, dimred="PCA")

################################################################################
# Graph-based clustering

print("Performing graph-based clustering on Harmony-corrected PCs.")

graph.clusters <- custom_clustering(reducedDim(clean,"harmony_PCA"), 25)
clean$graph <- paste(gsub("^","",LIB), unname(graph.clusters), sep="_")

################################################################################
# Write results

print(sprintf("Write filtered SingleCellExperiment object for L%s",LIB))

LIB_QC_FILE=sprintf("L%s_%scells_sce.RDS",LIB,ncol(clean))

saveRDS(clean,sprintf("%s/%s",LIB_QC_DIR,QC_FILE))
