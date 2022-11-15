################################################################################
################################################################################
# File name: 1b1__celltype_identification.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Four-step cell type assignment procedure
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
# Data loading

AGGR_QC_FILE="sce_rescaled_full_125libs_updated__harmony_run.RDS"

# load clean SingleCellExperiment object created in 1a2__quality_control__aggregated.R
rescaled=readRDS(sprintf("%s/data/%s",AGGR_QC_DIR,AGGR_QC_FILE))

################################################################################
# Data filtering

# add useful meta data columns
rescaled$Barcode_lib <- sprintf("%s-%s",rescaled$Barcode,rescaled$LIB)
rescaled$umap1 <- reducedDim(rescaled,"UMAP_harmony_run")[,1]
rescaled$umap2 <- reducedDim(rescaled,"UMAP_harmony_run")[,2]
rescaled$batch <- case_when(rescaled$run%in%as.character(1:14)~"Main",T~rescaled$run)
colnames(rescaled) <- rescaled$Barcode_lib

meta_data <- as.data.table(rescaled@colData)

# remove 1-cell 'clusters' induced after Harmony correction
nclusters=meta_data[,.N,by=cluster_seurat][N>5,max(cluster_seurat)]
clusters_1_cell=meta_data[cluster_seurat>nclusters,unique(cluster_seurat)]
rescaled <- rescaled[,ifelse(rescaled$cluster_seurat%in%clusters_1_cell,F,T)]

# remove cells from donors with < 500 cells in at least one condition
ind_lt_500cells=meta_data[,.N,by=.(IID,COND)][N<500,unique(IID)]
rescaled <- rescaled[,ifelse(rescaled$IID%chin%ind_lt_500cells,F,T)]

meta_data=meta_data[Barcode_lib%chin%rescaled$Barcode_lib]

# identify clusters with high proportion of flagged low-quality (5th percentile) barcodes

meta_data[,lib_QC5:=case_when(lib_QC5=="Pass"~T,lib_QC5=="Fail"~F)]
prop_loqc_barcodes=meta_data[,mean(!lib_QC5)*100,by=cluster_seurat]
setnames(prop_loqc_barcodes,"V1","prop_fail")

################################################################################
# Low-resolution metaclustering

de_design=meta_data[,.(cluster_seurat,Barcode_lib)][order(cluster_seurat),]

# find cluster-unconditional markers
de__clusters_ref_all <- lapply(X=unique(rescaled$cluster_seurat),FUN=function(x){
  # use cells from cluster x as query
  query=de_design[cluster_seurat==x,Barcode_lib]

  # if applicable, subsample to 10,000 observations
  if (length(query)>10000){
    query=query[sample(1:length(query),10000,F)]
  }

  # randomly sample the background (i.e., cells not in cluster x)
  background=de_design[cluster_seurat!=x,sample(Barcode_lib,10000,F)]

  # build the test data set and assign labels
  rescaled_sub <- rescaled[,c(query,background)]
  rescaled_sub$groups <- ifelse(rescaled_sub$Barcode_lib%chin%query,"query","background")

  # test differential expression (t-test for effect size, Wilcoxon for p-value)
  DE=findMarkers(rescaled_sub,groups=rescaled_sub$groups,test.type=c("wilcox"))@listData[["query"]]
  DE@listData$Symbol <- DE@rownames%>%{setNames(rowData(rescaled_sub)$Symbol,rowData(rescaled_sub)$ID)[.]}
  DE_FC=findMarkers(rescaled_sub,groups=rescaled_sub$groups)@listData[["query"]]

  # build output
  DE_FC=merge(
    data.table(ID=rownames(DE),as.data.frame(DE)),
    data.table(ID=rownames(DE_FC),as.data.frame(DE_FC)),
    by='ID',suffix=c('.rank','.norm')
  )

  DE_FC[,cluster_seurat:=x][,background:="ALL"]
  return(DE_FC)
})%>%rbindlist()

# define metaclusters based on CD14, CD3E, CD4, CD8A, CD19 and NCAM1 expression, by condition
T_NS_clusters=c(80,89,5,45,20,3,10,70,29,38,49,46,53,51,66,37,22,24,67,47,48,59,15,34,26,56,18,44,90)
T_S_clusters=c(0,1,11,42,25,28,64,2,17,82,54,12,16,4,13,68,52,69,57,9,6,36)
MONO_NS_clusters=c(23,74,71,41,55,14,32)
MONO_COV_clusters=c(7,33,39,60,65)
MONO_IAV_clusters=c(19,40,43,84)
DC_clusters=c(73,62,63)
PROG_clusters=c(76,78)
B_clusters=c(61,21,72,86,88,31,75,30,77,35,58,8,27,83,79)

meta_data[,metacluster:=case_when(
  cluster_seurat%in%c(T_NS_clusters,T_S_clusters)~"TNK",
  cluster_seurat%in%c(MONO_NS_clusters,MONO_COV_clusters,MONO_IAV_clusters,
    DC_clusters,PROG_clusters)~"MYELOID",
  cluster_seurat%in%c(B_clusters)~"B",
  T~"ALL"
)]
rescaled$metacluster=meta_data[,setNames(metacluster,Barcode_lib)][rescaled$Barcode_lib]

################################################################################
# High-resolution subclustering

# recluster cells within each metacluster and condition at higher resolution

metacluster_cluster=lapply(meta_data[,unique(metacluster)],function(m){
  round1=lapply(meta_data[,c("ALL",unique(COND))],function(c){

    if(c=="ALL"){
      cond=c("NS","COV","IAV")
    } else {
      cond=c
    }

    rescaled_sub=rescaled[,meta_data[metacluster==m&COND%chin%cond,Barcode_lib]]
    dec.pois.clean <- modelGeneVarByPoisson(rescaled_sub)
    rescaled_sub <- denoisePCA(rescaled_sub,technical=dec.pois.clean,name="PCA")
    reducedDim(rescaled_sub, "harmony_PCA_run") <- HarmonyMatrix(
      reducedDim(rescaled_sub,"PCA"),
      rescaled_sub$run,"run",do_pca=FALSE
    )
    rescaled_sub$cluster_seurat_sub <- custom_clustering(
      reducedDim(rescaled_sub,"harmony_PCA_run"),25,res=3
    )

    rescaled_sub <- runUMAP(rescaled_sub,dimred="harmony_PCA_run",name="harmony_UMAP_run")
    rescaled_sub$umap1_sub <- reducedDim(rescaled_sub,"harmony_UMAP_run")[,1]
    rescaled_sub$umap2_sub <- reducedDim(rescaled_sub,"harmony_UMAP_run")[,2]

    meta_data_sub=as.data.table(rescaled_sub@colData)
    meta_data_sub=meta_data_sub[,.(
      Barcode_lib,
      umap1=umap1_sub,umap2=umap2_sub,
      cluster=cluster_seurat,recluster=cluster_seurat_sub,
      metacluster=m,
      condition=c
    )]
    return(meta_data_sub)
  })
})

# find conditional markers given metacluster belonging

de_design=metacluster_cluster[,.(
  recluster=sprintf("%s_%s_%s",recluster,metacluster,condition),
  metacluster=sprintf("%s_%s",metacluster,condition),
  Barcode_lib)][order(metacluster),]

de__clusters_ref_metacluster <-
lapply(X=de_design[,unique(recluster)],FUN=function(x){
  # use cells from cluster x as query
  query=de_design[recluster==x,Barcode_lib]

  # if applicable, subsample to 10,000 observations
  if (length(query)>10000){
    query <- query[sample(1:length(query),10000,F)]
  }

  # randomly sample the background (i.e., cells not in cluster x but in the same meta cluster)
  mtclstr=de_design[recluster==x,unique(metacluster)]
  background=de_design[metacluster==mtclstr&recluster!=x,Barcode_lib]
  if (length(background)>10000){
    background <- background[sample(1:length(background),10000,F)]
  }

  # build the test data set and assign labels
  if (length(query>0)&length(background)>0){
  rescaled_sub <- rescaled[,c(query,background)]
  rescaled_sub$groups <- ifelse(rescaled_sub$Barcode_lib%in%query,"query","background")

  # test differential expression (t-test for effect size, Wilcoxon for p-value)
  DE=findMarkers(rescaled_sub,groups=rescaled_sub$groups,test.type=c("wilcox"))@listData[["query"]]
  DE@listData$Symbol <- DE@rownames %>% {setNames(rowData(rescaled_sub)$Symbol,rowData(rescaled_sub)$ID)[.]}
  DE_FC=findMarkers(rescaled_sub,groups=rescaled_sub$groups)@listData[["query"]]

  # build output
  DE_FC=merge(data.table(ID=rownames(DE),as.data.frame(DE)),data.table(ID=rownames(DE_FC),as.data.frame(DE_FC)),by='ID',suffix=c('.rank','.norm'))

  DE_FC[,recluster:=as.numeric(str_replace(x,"_.+$",""))][,background:=mtclstr]
  return(DE_FC)
  }
})%>%rbindlist()

# merge all unconditional and conditional markers
de__clusters <- rbind(de__clusters_ref_all,de__clusters_ref_metacluster)

fwrite(de__clusters,sprintf("%s/data/de__clusters.tsv.gz",DAT_DES_DIR),sep="\t")

meta_data[,updated_celltype:=make.names(as.character(celltype))]
meta_data[cluster_seurat==24,updated_celltype:='MAIT']
meta_data[cluster_seurat%in%c(40,60),updated_celltype:='MONO.CD14.INFECTED']
meta_data[cluster_seurat%in%c(83,88),updated_celltype:='B.INFECTED']
meta_data[cluster_seurat==4,updated_celltype:='T.CD8.EMRA']
meta_data[cluster_seurat==15,updated_celltype:='T.CD8.EMRA']
meta_data[cluster_seurat==9,updated_celltype:='NK.M.LIKE'] # KLRC2+ KLRC1- NK
meta_data[cluster_seurat==26,updated_celltype:='NK.M.LIKE']
meta_data[cluster_seurat==69,updated_celltype:='NK.INFECTED']
meta_data[cluster_seurat==64,updated_celltype:='T.N.INFECTED']
meta_data[cluster_seurat%in%c(59,48,57),updated_celltype:='MIX_NK.CD56dim_T.CD8.EMRA']
meta_data[cluster_seurat%in%c(90,52),updated_celltype:='NK.OTHER']
meta_data[cluster_seurat%in%c(42,45,47,48,80,51),updated_celltype:='MIX_T.CD4.N_T.CD8.N']
meta_data[cluster_seurat%in%c(72,75,86),updated_celltype:='B.OTHER']
meta_data[cluster_seurat%in%c(74),updated_celltype:='MONO.CD14.OTHER']
meta_data[cluster_seurat%in%c(46),updated_celltype:='T.CD4.E']
meta_data[cluster_seurat%in%c(82),updated_celltype:='ILC']
meta_data[cluster_seurat%in%c(49,68,67),updated_celltype:='T.OTHER']
meta_data[cluster_seurat%in%c(47,53,50,61,65,70,71,77,87),updated_celltype:='LowQC']

################################################################################
# Linear discriminant analysis

T_NS_clusters=c(80,89,5,45,20,3,10,70,29,38,49,46,53,51,66,37,22,24,67,47,48,59,15,34,26,56,18,44,90)
T_S_clusters=c(0,1,11,42,25,28,64,2,17,82,54,12,16,4,13,68,52,69,57,9,6,36)
MONO_NS_clusters=c(23,74,71,41,55,14,32)
MONO_COV_clusters=c(7,33,39,60,65)
MONO_IAV_clusters=c(19,40,43,84)
DC_clusters=c(73,62,63)
#PROG_clusters=c(76,78)
B_clusters=c(61,21,72,86,88,31,75,30,77,35,58,8,27,83,79,81,85)



meta_data$metacluster=case_when(
  meta_data$cluster_seurat%in%c(T_NS_clusters,T_S_clusters)~"TNK",
  meta_data$cluster_seurat%in%c(MONO_NS_clusters,MONO_COV_clusters,MONO_IAV_clusters,DC_clusters)~"MYELOID",
  meta_data$cluster_seurat%in%c(B_clusters)~"B",
  T~"ALL"
)

resolved_mixes=lapply(c("TNK","MYELOID","B"),function(mtc){
  resolved_mix_in_mtc=lapply(c("NS","COV","IAV"),function(cond){

    meta_cond=meta_data[metacluster==mtc&COND==cond,]

    mixes=meta_cond[grepl("^MIX",celltype_AB),unique(celltype_AB)]

    if (length(mixes)>0) {
    resolved_mix_in_cond=lapply(meta_cond[grepl("^MIX",celltype_AB),unique(celltype_AB)],function(MIX){
      A=str_split(MIX,"_",simplify=T)[,2]
      B=str_split(MIX,"_",simplify=T)[,3]
      # Whole data: cells from cell types A, B and MIX_A_B
      data=meta_cond[celltype_AB%in%c(A,B,MIX),]
      # Test data:  cells from MIX_A_B
      test_samples=data[celltype_AB==MIX,Barcode_lib]
      # Train data: cells from A and B
      subset=5000
      A_samples=data[celltype_AB==A,Barcode_lib]
      B_samples=data[celltype_AB==B,Barcode_lib]
      A_prop=round(length(A_samples)/length(c(A_samples,B_samples)),2)
      B_prop=round(length(B_samples)/length(c(A_samples,B_samples)),2)
      A_num=round(subset*A_prop)
      B_num=round(subset*B_prop)
      if (length(c(A_samples,B_samples))>subset){
        A_samples=data[celltype_AB==A,][sample(1:.N,A_num,replace=F),Barcode_lib]
        B_samples=data[celltype_AB==B,][sample(1:.N,B_num,replace=F),Barcode_lib]
      } else {
        A_samples=data[celltype_AB==A,Barcode_lib]
        B_samples=data[celltype_AB==B,Barcode_lib]
      }

      train_samples=c(A_samples,B_samples)

      # Build train data: keep 10% most variable genes
      train_data=logcounts(sce[,train_samples])
      features_to_keep=setNames(rowVars(train_data),row.names(train_data))%>%sort(decreasing=T)%>%{names(.[1:round(12667*0.1)])}
      train_data=as.data.table(t(scale(train_data[features_to_keep,])))
      rownames(train_data)=train_samples
      train_data$celltype_AB=data[,setNames(celltype_AB,Barcode_lib)][train_samples]

      # Fit the model
      model=lda(celltype_AB~., data=train_data)
      # Buil test data
      test_data=as.data.table(t(scale(logcounts(sce[features_to_keep,test_samples]))))

      # Predict cell types in mixed cluster
      predictions=predict(model,test_data)
      delta_posterior=abs(predictions$posterior[,1]-predictions$posterior[,2])
      resolved_mix=data.table(Barcode_lib=test_samples,celltype_resolved=predictions$class,delta_posterior)
      resolved_mix[,mix:=MIX]
      resolved_mix[,metacluster:="TNK"]
      return(resolved_mix)
      })%>%{do.call("rbind",.)}
      resolved_mix_in_cond[,condition:=cond]
     return(resolved_mix_in_cond)
     }
  })%>%{do.call("rbind",.)}
  resolved_mix_in_mtc[,metacluster:=mtc]
 return(resolved_mix_in_mtc)
})%>%{do.call("rbind",.)}

resolved_mixes$celltype_resolved_confident=ifelse(
  resolved_mixes$delta_posterior>0.80,
  as.character(resolved_mixes$celltype_resolved),
  "MixUnresolved"
)

meta_data$celltype_AB_resolved=resolved_mixes[,setNames(as.character(celltype_resolved_confident),Barcode_lib)][meta_data$Barcode_lib]
meta_data$celltype_AB_resolved=ifelse(
  is.na(meta_data$celltype_AB_resolved),
  meta_data$celltype_AB,
  meta_data$celltype_AB_resolved
)

################################################################################
# Merge results and write

meta_data=merge(
  meta_data,
  metacluster_cluster[condition!="ALL",.(Barcode_lib,lores_cluster=cluster,hires_cluster=recluster)],
  by="Barcode_lib",all.x=T
)

meta_data$celltype_AB_resolved=ifelse(
  meta_data$lores_cluster%in%c(40,60),
  "MONO.CD14.INFECTED",
  meta_data$celltype_AB_resolved
)

meta_data$celltype_AB_resolved=ifelse(
  meta_data$celltype_AB_resolved=="mDC",
  "MONO.CD14",
  meta_data$celltype_AB_resolved
)

meta_data$metacluster=ifelse(
  meta_data$celltype_AB_resolved=="Plasmablast",
  "B",
  meta_data$metacluster
)

fwrite(meta_clean,sprintf("%s/data/sce_clean_metadata.tsv",AGGR_QC_DIR),sep="\t")
