################################################################################
################################################################################
# File name: 1e3__scores.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Compute ISG and inflammatory scores per lineage, donor and condition
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
# ISG and inflammatory score by IID, lineage and condition

# reshape batch-adjusted counts and load onto Seurat object
Expr_matrix=dcast(Expr,ID~IID+celltype+state,value.var="logCPM")%>%as.matrix(rownames="ID")
seurat=CreateSeuratObject(counts=Expr_matrix)
seurat$IID=str_split(colnames(seurat),"_",simplify=T)[,1]
seurat$celltype=str_split(colnames(seurat),"_",simplify=T)[,2]
seurat$condition=str_split(colnames(seurat),"_",simplify=T)[,3]

# extract feature meta data
feature_data=Expr[,.(ID,Symbol)]%>%unique

REF_DIR="1__transcriptome_processing/data"

# load hallmark interferon-stimulated pathway genes
hallmark_ifng_response_genes=fread(sprintf("%s/hallmark_ifng_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
hallmark_ifna_response_genes=fread(sprintf("%s/hallmark_ifna_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
hallmark_ifn_response_genes=unique(c(hallmark_ifng_response_genes,hallmark_ifna_response_genes))
names(hallmark_ifn_response_genes)=unique(DE_COND[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_ifn_response_genes]

# load hallmark inflammatory pathway genes
hallmark_inflammatory_response_genes=fread(sprintf("%s/hallmark_inflammatory_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
names(hallmark_inflammatory_response_genes)=unique(DE_COND[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_inflammatory_response_genes]

hallmark_inflammatory_exclusive=which(!hallmark_inflammatory_response_genes%in%hallmark_ifn_response_genes)%>%hallmark_inflammatory_response_genes[.]
hallmark_ifn_exclusive=which(!hallmark_ifn_response_genes%in%hallmark_inflammatory_response_genes)%>%hallmark_ifn_response_genes[.]

# run AddModuleScores on Seurat object
seurat=AddModuleScore(seurat,
  features=list(
    hallmark_ifna_response_score=names(hallmark_ifna_response_genes),
    hallmark_ifng_response_score=names(hallmark_ifng_response_genes),
    hallmark_isg_score=names(hallmark_ifn_exclusive),
    hallmark_inflammatory_score=names(hallmark_inflammatory_exclusive)
  ),
  name=c("hallmark_ifna_response_score","hallmark_ifng_response_score","hallmark_isg_score","hallmark_inflammatory_score","go_ifna_response_score")
)
colnames(seurat@meta.data) %<>% str_replace("[0-9]+$","")

as.data.table(seurat@meta.data)[,-c("orig.ident","nCount_RNA","nFeature_RNA")]%>%
fwrite(.,"1__transcriptome_processing/data/scores__by__IID_lineage_condition.tsv",sep="\t")

################################################################################
# ISG and inflammatory score by IID, lineage and condition

Expr_matrix=dcast(Expr,ID~IID+state,value.var="logCPM",fill=0)%>%as.matrix(rownames="ID")

seurat=CreateSeuratObject(counts=Expr_matrix)
seurat$IID=str_split(colnames(seurat),"_",simplify=T)[,1]
seurat$celltype=str_split(colnames(seurat),"_",simplify=T)[,2]
seurat$condition=str_split(colnames(seurat),"_",simplify=T)[,3]

# run AddModuleScores on Seurat object
seurat=AddModuleScore(seurat,
  features=list(
    hallmark_ifna_response_score=names(hallmark_ifna_response_genes),
    hallmark_ifng_response_score=names(hallmark_ifng_response_genes),
    hallmark_isg_score=names(hallmark_ifn_exclusive),
    hallmark_inflammatory_score=names(hallmark_inflammatory_exclusive)
  ),
  name=c("hallmark_ifna_response_score","hallmark_ifng_response_score","hallmark_isg_score","hallmark_inflammatory_score","go_ifna_response_score")
)
colnames(seurat@meta.data) %<>% str_replace("[0-9]+$","")

as.data.table(seurat@meta.data)[,-c("orig.ident","nCount_RNA","nFeature_RNA")]%>%
fwrite(.,"1__transcriptome_processing/data/scores__by__IID_condition.tsv",sep="\t")
