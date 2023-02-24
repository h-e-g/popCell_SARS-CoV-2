#!/bin/bash
################################################################################
################################################################################
# File name: 3a9__launcher_part_a5_a6_a7.sh
# Author: Y.A., M.R., M.ON.
################################################################################
###############################################################################
# Step: aggregate celltype/lineage eQTLs, then retrieve and combine eQTL sumstats from all chromosomes/celltypes/conditions
# the script will do aggregation for both eQTL and reQTLs
# Batch launcher script
################################################################################
################################################################################

################################################################################
#### Example calls
#./3a9__launcher_part_a54_a6_a7.sh

LOGS="./LOGS"
JOBFILE='./LOGS/launched_jobs_3a9.txt'

MISC_DIR='./MISC/'
EQTL_SCRIPT_DIR="./3__eQTL_mapping/3a_eQTL_mapping"
EQTL_DIR="3__eQTL_mapping/sumStats"

DIST=100000
RUNID=RUN1

################################# ################################# #################################
#################################         step 1: merging           #################################
################################# ################################# #################################

OPTIONS_MERGE="--dist ${DIST}"
SBATCH_MERGE_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${LOGS}/logfile_eQTL_merge_celltype_lineage_%A_%a.log -J eQTL_merge"

################################# merge lineage and celltype eQTLs
RUN_NAME=lineage_condition___CellPropLineage_SVs_RUN1
RUN_NAME_CELLTYPE=celltype_condition___CellPropLineage_SVs_RUN1

JOBID_MERGE=`sbatch ${SBATCH_MERGE_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a5__merge_eQTL_celltype_and_lineage.R ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE} --run_id ${RUNID}_eQTL`
echo "${JOBID_MERGE} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}" >> ${JOBFILE}
echo ""

################################# merge lineage and celltype reQTLs
RUN_NAME=lineage_condition_logFC__CellPropLineage_SVs_RUN1
RUN_NAME_CELLTYPE=celltype_condition_logFC__CellPropLineage_SVs_RUN1

JOBID_MERGE_LOGFC=`sbatch ${SBATCH_MERGE_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a5__merge_eQTL_celltype_and_lineage.R ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE} --run_id ${RUNID}_reQTL`
echo "${JOBID_MERGE_LOGFC} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}" >> ${JOBFILE}
echo ""

################################# ################################# #################################
################################# step 2: add stats per chromosome  #################################
################################# ################################# #################################

# followed by

################################# ################################# #################################
################################# step 3: run aggregate on all chr  #################################
################################# ################################# #################################

OPTIONS="--dist ${DIST}"
SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G -o ${LOGS}/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
SBATCH_AGG_ARGS="--parsable --mem=100G -o ${LOGS}/logfile_eQTL_AGG_allCT_%A_%a.log -J eQTL_AGG"

############@ Annotate lineage & cell type level eQTLs

##############
# step 2 : add lineage expr stats
RUN_NAME=lineage_condition___CellPropLineage_SVs_RUN1
 JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} --dependency=afterok:${JOBID_MERGE} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a6__retrieve_eQTLSumStats_allCond_1CHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_eQTL --chr`
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id '${RUNID}_eQTL' --chr 1-22" >> ${JOBFILE}
 echo ""

# step 3: aggregate lineage expr stats
  JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a7__aggregate_eQTLSumStats_allCond_allCHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_eQTL`
  echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id '${RUNID}_eQTL'" >> ${JOBFILE}
  echo ""


##############
# step 2 :add celltype expr stats
 RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_RUN1
 JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} --dependency=afterok:${JOBID_MERGE} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a6__retrieve_eQTLSumStats_allCond_1CHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_eQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr`
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id 'RUN1_eQTL' --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
 echo ""

# step3: aggregate celltype expr stats
RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_RUN1
JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a7__aggregate_eQTLSumStats_allCond_allCHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_eQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype`
echo "${JOBID_AGG} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_eQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
echo ""


 ############@ Annotate lineage & cell type level reQTLs
 RUN_NAME=lineage_condition_logFC__CellPropLineage_SVs_RUN1
 #step 2: add lineage logFC stats
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} --dependency=afterok:${JOBID_MERGE_LOGFC} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a6__retrieve_eQTLSumStats_allCond_1CHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL --chr 1-22" >> ${JOBFILE}
echo ""

# step3: aggregate lineage logFC stats
JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a7__aggregate_eQTLSumStats_allCond_allCHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL`
 echo "${JOBID_AGG} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL" >> ${JOBFILE}
 echo ""


##############
# step2: add celltype logFC stats
RUN_NAME_ANNOT=celltype_condition_logFC__CellPropLineage_SVs_RUN1
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} --dependency=afterok:${JOBID_MERGE_LOGFC} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a6__retrieve_eQTLSumStats_allCond_1CHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id 'RUN1_reQTL' --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
echo ""


# step3: celltype logFC stats
JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a7__aggregate_eQTLSumStats_allCond_allCHR.R ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype`
echo "${JOBID_AGG} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_id ${RUNID}_reQTL --run_name_annotate ${RUN_NAME_ANNOT} --name celltype" >> ${JOBFILE}
echo ""

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
