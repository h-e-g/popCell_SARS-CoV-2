######## objective of the script:
# for each level of CELLTYPE and STATE, compute count and TPM for every gene and individual
# for each STATE, compute the percentage of various CELLTYPEs for all individuals
# 		(allows assessing the impact of cellular interactions: eg. how % of Monocytes affects T cell response)
# for each STATE and CELLTYPE, compute the percentage of various subsets (seurat clusters) for all individuals
# 		(allows assessing the impact of cellular heterogeneity:  how % CD8+ EM cells affect the overall T cell response)

################################################################################
# START OBSOLETE: we now work with a unique 21-level cell type factor

###### requirements:
### define CELLTYPE and STATE
# ideally CELLTYPE and CELLSTATE, should not contain any '_'
# eg. CELLTYPE='celltype.19level'
#     STATE='activation.state'
#
# or CELLTYPE='celltype.8level'
#     STATE='condition'
#
# STOP OBSOLETE
################################################################################


### define path to an sceObject (SCE_OBJECT) containing :
# only 6 hour time point,
# only NS,COV and IAV conditions
# only individuals with >400 cells in all 3 conditions (and for which activation is reproducible)
# undefined/Dying cells should be have been removed already

### we expect the following variables to be defined in the sceobject
# cluster_seurat, condition, celltype, condition

### set NLIBS to the number of libraries in the current iteration (for naming purposes)

### set OUT_DIR to the desired outpout directory ( evo_immuno_pop/single_cell/project/pop_eQTL/data by default)

###### start here
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

# default parameter
SCE_OBJECT='sce_clean.RDS'
NLIBS=125 # number of libraries in the dataset (for naming purpose only)
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition' # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
OUT_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences",EVO_IMMUNO_POP_ZEUS) # output directory (defaults to single_cell/project/pop_eQTL/data )

# update parameter values based on provided arguments
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--sce' | cmd[i]=='-s' ){SCE_OBJECT = cmd[i+1]} # path to the sce object to take as input
	if (cmd[i]=='--count' | cmd[i]=='-c' ){COUNT_DATA_FILE = cmd[i+1]} # path to count & CPM data file withs one line per LIB, IID, CELLTYPE, STATE, and gene.
	if (cmd[i]=='--state' | cmd[i]=='-a' ){STATE = cmd[i+1]} # state variable to use (cellular activation state or experimental condition). Will be used for naming of output files
  if (cmd[i]=='--nlibs' | cmd[i]=='-n' ){NLIBS = cmd[i+1]} # number of libraries in the dataset (for naming purpose only)
  if (cmd[i]=='--outdir' | cmd[i]=='-o' ){OUT_DIR = cmd[i+1]} # output directory (defaut: single_cell/project/pop_eQTL/data/2_population_differences)
}


sce.Object=readRDS(sprintf("%s/%s",DATA_DIR,SCE_OBJECT))
sce.Object$condition <- sce.Object$COND

# Add lineage information

meta_data_cells=as.data.table(sce.Object@colData)
setnames(meta_data_cells,'cluster_seurat','cluster')
setnames(meta_data_cells,STATE,'state')

################################################################################
##### Aggregate data into pseudobulk by library, IID, cell type and condition ##
################################################################################

ColDATA=sce.Object@colData[,c(STATE,"IID","LIB")]

print(str(ColDATA))

summed <- scuttle::aggregateAcrossCells(sce.Object, id=ColDATA)
saveRDS(summed,file=sprintf('%s/sce_pseudoBulk_full_%slibs_per_%s_IID_LIB.RDS',OUT_DIR,NLIBS,STATE))

Counts=melt(data.table(ID=rownames(counts(summed)),counts(summed)),value.name='count',variable.name='group_id')
sample_annot=as.data.table(colData(summed))
setnames(sample_annot,STATE,'state')
sample_annot=sample_annot[,.(IID, LIB, ncells, state)]
sample_annot[,'group_id':=paste('V',1:.N,sep='')]
Counts=merge(Counts,sample_annot,by='group_id')
gene_annot=as.data.table(rowData(summed))[,.(ID,Symbol)]
Counts=merge(Counts,gene_annot,by='ID')
Counts[,CPM:=count/sum(count)*1e6,by=group_id]
Counts[,group_id:=NULL]
fwrite(Counts,file=sprintf('%s/Counts_and_CPM__pseudoBulk_full_%slibs__per_%s_IID_and_LIB.tsv.gz',OUT_DIR, NLIBS, STATE))
