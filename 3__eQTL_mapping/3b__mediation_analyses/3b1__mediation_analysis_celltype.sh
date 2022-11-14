#!/bin/bash
################################################################################
################################################################################
# File name: 3b1__mediation_analysis_celltype.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use mediation analysis to quantify var. explained by cell proportions
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=3b1
#SBATCH --output=${LOG_DIR}/3b1__%J.log
#SBATCH --mem=100G
################################################################################

################################################################################
# Setup

# load modules
module load R/4.1.0

################################################################################
# Command history

# fraction of gene expression variance in each condition

# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage B --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE

# sbatch ./mediation_analysis_celltype.sh --state NS --lineage MONO --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage B --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage T.CD4 --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage T.CD8 --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage NK --comparison ASH_AFB --runid 220409

# sbatch ./mediation_analysis_celltype.sh --state NS --lineage MONO --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage B --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage T.CD4 --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage T.CD8 --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_celltype.sh --state NS --lineage NK --comparison ASH_EUB --runid 220409

# fraction of gene expression log-fold change variance in each stimulated condition

# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state IAV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR

# sbatch ./mediation_analysis_celltype.sh --state COV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state COV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state COV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state COV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_celltype.sh --state COV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR

Rscript ./3b1__mediation_analysis_celltype.R $*

exit
