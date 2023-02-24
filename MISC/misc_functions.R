# aggregate SingleCellExperiment objects
CombineSCE=function(SCE_list){
  for(i in 1:length(SCE_list)){
    reducedDims(SCE_list[[i]])=NULL
    metadata(SCE_list[[i]]) <- list()
    SCE_list[[i]]$sizeFactor <- NULL
    names(rowData(SCE_list[[i]])) <- c("ID","Symbol","Type")
    SCE_list[[i]]$Barcode=paste(SCE_list[[i]]$Barcode,substr(names(SCE_list)[i],2,100),sep='-')
    colnames(SCE_list[[i]]) <- SCE_list[[i]]$Barcode
  }
  do.call(cbind,SCE_list)
}

# graph-based clustering Seurat's way
custom_clustering <- function(pcs,k.param,res=0.8){
nn.ranked <- RANN::nn2(data=pcs,k=k.param,searchtype="standard",eps=0)
nn.ranked <- nn.ranked$nn.idx
j <- as.numeric(t(nn.ranked))
i <- ((1:length(j))-1) %/% k.param + 1
nn.matrix <- as(
  Matrix::sparseMatrix(
    i=i,j=j,x=1,
    dims=c(nrow(pcs),nrow(pcs))
  ),
  Class="Graph")
rownames(nn.matrix) <- rownames(pcs)
colnames(nn.matrix) <- rownames(pcs)
neighbor.graphs <- list(nn=nn.matrix)
# Compute shared nearest neighbors
prune.SNN <- 1/15
ComputeSNN <- function(nn_ranked, prune) {
    .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn_ranked, prune)
}
snn.matrix <- ComputeSNN(
  nn_ranked = nn.ranked,
  prune = prune.SNN
)
rownames(snn.matrix) <- rownames(pcs)
colnames(snn.matrix) <- rownames(pcs)
snn.matrix <- as.Graph(x = snn.matrix)
neighbor.graphs[["snn"]] <- snn.matrix

RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
    .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat', SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
}
edge.file.name <- NULL
edge_file <- ''
k25_snn_louvain <- RunModularityClusteringCpp(
  SNN=neighbor.graphs[["snn"]],
  modularityFunction=1,
  resolution=res,
  algorithm=1,
  nRandomStarts=10,
  nIterations=10,
  randomSeed=0,
  printOutput=T,
  edgefilename=''
)
names(k25_snn_louvain) <- colnames(neighbor.graphs$snn)
return(k25_snn_louvain)
}

# rank transformation
rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
}

# not in
`%nin%`=Negate(`%in%`)

# adjust a variable given covariates
adjust_on_x=function(y,matx,matz){
    # y : variable to adjust (a vector)
  # mat_x: a vector, matrix or data.frame containing variables to adjust on
  # mat_z:  a vector, matrix or data.frame containing variables to be included in the model and not adjusted on
  if(is.data.frame(matx)){matx=as.matrix(matx)}
  if(is.vector(matx)){matx=matrix(matx,ncol=1)}
  if(is.data.frame(matz)){matz=as.matrix(matz)}
  if(is.vector(matz)){matz=matrix(matz,ncol=1)}

    y_z=lm(y~.,data=as.data.frame(matz))$res
    matx_z=apply(matx,2,function(x){lm(x~.,data=as.data.frame(matz))$res})
    beta_x_z=lm(y_z~.,data=as.data.frame(matx_z))$coef[-1]
    y_adj=y-scale(matx,T,F)%*%matrix(beta_x_z,ncol=1)
    as.vector(y_adj)
}
