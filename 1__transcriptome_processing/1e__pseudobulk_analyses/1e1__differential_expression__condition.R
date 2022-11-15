################################################################################
################################################################################
# File name: 1e1__differential expression__condition.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Condition differential expression tests
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
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition'

# load batch-adjusted counts computed in 1c2__pseudobulk_batch_correction.R
Expr_ct=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))
CELLTYPE='lineage'
Expr_ln=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))

################################################################################
# Differential expression between conditions: lineage-level

# reshape batch-adjusted pseudobulk data
Expr_test=dcast(Expr_ln,IID+celltype+Symbol+ID~state,value.var="logCPM")

# difference in means, Wilcoxon's p-value by gene, lineage and condition
DE_COND=Expr_test[,.(logFC_COV=mean(COV)-mean(NS),
                  pval_COV=wilcox.test(COV-NS,data=.SD)$p.value,
                  logFC_IAV=mean(IAV)-mean(NS),
                  pval_IAV=wilcox.test(IAV-NS,data=.SD)$p.value)
                  ,by=.(ID,Symbol,celltype)]

DE_COND[,FDR_COV:=p.adjust(pval_COV,"fdr")]
DE_COND[,FDR_IAV:=p.adjust(pval_IAV,"fdr")]

# keep only significantly differentially expressed genes
logfcthr=0.5
DE_COND[grepl('ENSG', ID) & ((FDR_COV<0.01 & logFC_COV>logfcthr)|(FDR_IAV<0.01 & logFC_IAV>logfcthr)),length(unique(Symbol))]

fwrite(DE_COND,"1__transcriptome_processing/data/DE_condition_by_lineage.tsv.gz",sep="\t")

################################################################################
# Differential expression between conditions: cell type-level

# reshape batch-adjusted pseudobulk data
Expr_test=dcast(Expr_ct,IID+celltype+Symbol+ID~state,value.var="logCPM")

# difference in means, Wilcoxon's p-value by gene, cell type and condition
DE_COND=Expr_test[,.(logFC_COV=mean(COV)-mean(NS),
                  pval_COV=wilcox.test(COV-NS,data=.SD)$p.value,
                  logFC_IAV=mean(IAV)-mean(NS),
                  pval_IAV=wilcox.test(IAV-NS,data=.SD)$p.value)
                  ,by=.(ID,Symbol,celltype)]

DE_COND[,FDR_COV:=p.adjust(pval_COV,"fdr")]
DE_COND[,FDR_IAV:=p.adjust(pval_IAV,"fdr")]

# keep only significantly differentially expressed genes
logfcthr=0.5
DE_COND[grepl('ENSG', ID) & ((FDR_COV<0.01 & logFC_COV>logfcthr)|(FDR_IAV<0.01 & logFC_IAV>logfcthr)),length(unique(Symbol))]

fwrite(DE_COND,"1__transcriptome_processing/data/DE_condition_by_celltype.tsv.gz",sep="\t")
