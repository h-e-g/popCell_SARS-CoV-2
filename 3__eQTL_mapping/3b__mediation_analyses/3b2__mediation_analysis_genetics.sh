#!/bin/bash
################################################################################
################################################################################
# File name: 3b1__mediation_analysis_genetics.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use mediation analysis to quantify var. explained by best eQTL
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"

################################################################################
# SLURM options
#SBATCH --job-name=3b2
#SBATCH --output=${LOG_DIR}/3b1__%J.log
#SBATCH --mem=100G
################################################################################

################################################################################
# Setup

# load modules
module load R/4.1.0
module load samtools

################################################################################
# Command history

# fraction of gene expression variance in each condition

# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage B --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff FALSE --exp_resp DE

# sbatch ./mediation_analysis_genetics.sh --state NS --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state NS --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state NS --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state NS --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state NS --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
#
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
#
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DE

# sbatch ./mediation_analysis_genetics.sh --state COV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
#
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage MONO --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage B --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD4 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD8 --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage NK --comparison EUB_AFB --runid 220617 --popdiff TRUE --exp_resp DR

# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage MONO --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage B --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD4 --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage T.CD8 --comparison ASH_AFB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state IAV --lineage NK --comparison ASH_AFB --runid 220409

# sbatch ./mediation_analysis_genetics.sh --state COV --lineage MONO --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage B --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD4 --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage T.CD8 --comparison ASH_EUB --runid 220409
# sbatch ./mediation_analysis_genetics.sh --state COV --lineage NK --comparison ASH_EUB --runid 220409

Rscript ./3b2__mediation_analysis_genetics.R $*

exit
