
################################################################################
################################################################################
# File name: 3a10_extract_eQTL_genotypes.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: create a tsv file with genotype data from all eQTLs/reQTLs aross 222 individuals
# This aims to allow faster loeading of eQTL genotypes down the line (eg. when plotting eQTLs)
# Effector script
#################################################################
#################################################################

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
EQTL_DIR = "../../../3__eQTL_mapping/SumStats"

source(sprintf("%s/load_eQTLs.R",MISC_DIR))
all_eQTL_snps=unique(c(eQTL_Signif_both$snps,reQTL_Signif_both$snps))


source(sprintf("%s/querySNPs.R",MISC_DIR))
SNP_info=getMap(annotate=TRUE)

tic('loading eQTL genotypes')
eQTL_genotypes=getSNP(all_eQTL_snps,vector=FALSE, Map=SNP_info)
eQTL_genotypes[,IID:=as.character(IID)]
toc()

keptIID=fread(sprintf('%s/data/IID_individual_to_use.tsv',DAT_EQTL_DIR),header=F)[,V1]
eQTL_genotypes=eQTL_genotypes[ID%chin%all_eQTL_snps & IID%chin%keptIID,]
fwrite(eQTL_genotypes,file=sprintf('%s/All_eQTL_and_reQTL_genotypes.tsv.gz',EQTL_DIR),sep='\t')
