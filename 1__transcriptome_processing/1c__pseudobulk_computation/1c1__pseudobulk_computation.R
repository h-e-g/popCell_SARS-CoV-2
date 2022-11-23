################################################################################
################################################################################
# File name: 1c1__pseudobulk_computation.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Pseudobulk computation from filtered SingleCellExperiment object
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
SCE_OBJECT='sce_clean.RDS'
NLIBS=125
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition'
OUT_DIR = "1__transcriptome_processing" # output directory

# update parameter values based on arguments provided in 1c__pseudobulk_computation.sh
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
        if (cmd[i]=='--sce' | cmd[i]=='-s' ){SCE_OBJECT = cmd[i+1]} # path to the sce object to take as input
        if (cmd[i]=='--count' | cmd[i]=='-c' ){COUNT_DATA_FILE = cmd[i+1]} # path to count & CPM data file withs one line per LIB, IID, CELLTYPE, STATE, and gene.
  if (cmd[i]=='--celltype' | cmd[i]=='-t' ){CELLTYPE = cmd[i+1]} # number of Surrogate variables to include (-1: determinaed automatically, 0: none, n>0 : only use the first n SVs )
        if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
  if (cmd[i]=='--nlibs' | cmd[i]=='-n' ){NLIBS = cmd[i+1]} # number of libraries in the dataset (for naming purpose only)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output directory (defaut: single_cell/project/pop_eQTL/data/2_population_differences)
}

################################################################################
# Data loading

AGGR_QC_FILE="sce_rescaled_full_125libs_updated__harmony_run.RDS"

# load clean SingleCellExperiment object created in 1a2__quality_control__aggregated.R
sce.Object=readRDS(sprintf("%s/data/%s",AGGR_QC_DIR,AGGR_QC_FILE))

# load meta data created in 1b1__celltype_identification.R
meta_clean=fread(sprintf("%s/data/sce_clean_metadata.tsv",AGGR_QC_DIR))

sce.Object@colData <- meta_clean
sce.Object$condition=sce.Object$COND

# add lineage information
if (CELLTYPE!="celltype") {
  celltype_lineage=fread("0__barcode_processing/data/lineages_celltype.tsv")
  sce.Object$lineage=celltype_lineage[,setNames(lineage,celltype)][sce.Object$celltype]
  sce.Object$celltype=NULL
}

meta_data_cells=as.data.table(sce.Object@colData)
setnames(meta_data_cells,'cluster_seurat','cluster')
setnames(meta_data_cells,CELLTYPE,'celltype')
setnames(meta_data_cells,STATE,'state')


dir.create(sprintf('%s/Cellcounts/',OUT_DIR))
dir.create(sprintf('%s/Cellcounts/Pct_cluster_by_%s_%s/',OUT_DIR,CELLTYPE,STATE))

# calculate cell type frequencies (%) by donor and condition
Cell_count=meta_data_cells[,.N,by=.(celltype,state,IID)]
Cell_count[,Ntot:=sum(N),by=.(celltype,state)]
Cell_count[,Pct_total:=N/sum(N),by=.(IID,state)]
Cell_count=dcast(Cell_count,IID+state~ celltype,value.var='Pct_total',fill=0)
fwrite(Cell_count,file=sprintf('%s/Cellcounts/Pct_%s_by_IID_%s.tsv.gz',OUT_DIR,CELLTYPE,STATE),sep='\t')

# calculate cluster frequencies (%) by cell type, donor and condition
Cell_count_cluster=meta_data_cells[,.N,by=.(cluster,celltype,state,IID)]
Cell_count_cluster[,Ntot:=sum(N),by=.(celltype,state)]
Cell_count_cluster[,Pct_broad:=N/sum(N),by=.(IID,celltype,state)]

groups=Cell_count_cluster[,.N,by=.(celltype,state)]
for (i in 1:groups[,.N]){
  myCELLTYPE=groups[i,celltype]
  mySTATE=groups[i,state]
  Pct_subsets_cluster=dcast(Cell_count_cluster[celltype==myCELLTYPE & state==mySTATE,],IID~paste('Clust',cluster,sep=''),value.var='Pct_broad',fill=0)
  fwrite(Pct_subsets_cluster,file=sprintf('%s/Cellcounts/Pct_cluster_by_%s_%s/Pct_cluster_by_IID_in_%s_%s.tsv.gz', OUT_DIR, CELLTYPE, STATE, myCELLTYPE, mySTATE),sep='\t')

################################################################################
# Compute pseudobulks by library, donor, cell type and condition

ColDATA=sce.Object@colData[,c(CELLTYPE,STATE,"IID","LIB")]

print(str(ColDATA))

summed <- scuttle::aggregateAcrossCells(sce.Object, id=ColDATA)
saveRDS(summed,file=sprintf('%s/data/sce_pseudoBulk_full_%slibs_per_%s_%s_IID_LIB.RDS',OUT_DIR,NLIBS,CELLTYPE,STATE))

Counts=melt(data.table(ID=rownames(counts(summed)),counts(summed)),value.name='count',variable.name='group_id')
sample_annot=as.data.table(colData(summed))
setnames(sample_annot,CELLTYPE,'celltype')
setnames(sample_annot,STATE,'state')
sample_annot=sample_annot[,.(IID, LIB, ncells, celltype, state)]
sample_annot[,'group_id':=paste('V',1:.N,sep='')]
Counts=merge(Counts,sample_annot,by='group_id')
gene_annot=as.data.table(rowData(summed))[,.(ID,Symbol)]
Counts=merge(Counts,gene_annot,by='ID')
Counts[,CPM:=count/sum(count)*1e6,by=group_id]
Counts[,group_id:=NULL]

fwrite(Counts,file=sprintf('%s/data/pseudobulk_%slibs__per_%s_%s_IID_and_LIB.tsv.gz',OUT_DIR, NLIBS, CELLTYPE, STATE))
