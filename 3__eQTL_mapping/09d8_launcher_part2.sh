#!/bin/bash

#./09d8_launcher_part2.sh


JOBFILE=zz_joblist_25052022.txt
SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
EQTL_SCRIPT_DIR="09_eQTLmapping"
EQTL_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/3_eQTL_mapping"

DIST=100000
USE_JOINT=FALSE
USE_CELLTYPE=TRUE

CORRECT_CELLTYPE=FALSE
CONDITIONAL=FALSE



################################# ################################# #################################
#################################         step 1: merging           #################################
################################# ################################# #################################

################################# merge lineage and celltype eQTLs

RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
RUN_NAME_CELLTYPE=celltype_condition___CellPropLineage_SVs_220409

OPTIONS_MERGE="--dist ${DIST} --add_joint_mapping ${USE_JOINT}"
SBATCH_MERGE_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_merge_celltype_lineage_%A_%a.log -J eQTL_merge"
JOBID_MERGE=`sbatch ${SBATCH_MERGE_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d5_merge_eQTL_celltype_and_lineage.R ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}`
echo "${JOBID_MERGE} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}"
echo "${JOBID_MERGE} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}" >> ${JOBFILE}
echo ""

################################# merge lineage and celltype reQTLs
RUN_NAME=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
RUN_NAME_CELLTYPE=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409

OPTIONS_MERGE="--dist ${DIST} --add_joint_mapping ${USE_JOINT}"
SBATCH_MERGE_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_merge_celltype_lineage_%A_%a.log -J eQTL_merge"
JOBID_MERGE_LOGFC=`sbatch ${SBATCH_MERGE_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d5_merge_eQTL_celltype_and_lineage.R ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}`
echo "${JOBID_MERGE} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}"
echo "${JOBID_MERGE} submitted ${OPTIONS_MERGE} --run_name ${RUN_NAME} --run_name_celltype ${RUN_NAME_CELLTYPE}" >> ${JOBFILE}
echo ""

################################# ################################# #################################
################################# step 2: add stats per chromosome  #################################
################################# ################################# #################################

# followed by

################################# ################################# #################################
################################# step 3: run mash on all per chr   #################################
################################# ################################# #################################

OPTIONS="--dist ${DIST} --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR --dependency=afterok:${JOBID_MERGE}"

############@ Annotate lineage & cell types expression

##############
# step 2 :add lineage expr stats
RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
 JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --chr`
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22"
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22" >> ${JOBFILE}
 echo ""

# step 3: mash lineage expr stats
SBATCH_AGG_MASH_ARGS="--parsable --mem=100G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_mash_%A_%a.log -J eQTL_AGG_mash"
  JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME}`
  echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22"
  echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22" >> ${JOBFILE}
  echo ""



##############
# step 2 :add celltype expr stats
 RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
 JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr`
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22"
 echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
 echo ""

# step3: mash celltype expr stats
RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
echo ""



##############
#step2: add lineage logFC stats
RUN_NAME_ANNOT=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22" >> ${JOBFILE}
echo ""

# step3: mash lineage logFC stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC" >> ${JOBFILE}
echo ""


##############
#step2: add celltype logFC stats
RUN_NAME_ANNOT=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC_celltype --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22" >> ${JOBFILE}
echo ""

# step3: mash celltype logFC stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC_celltype`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC " >> ${JOBFILE}
echo ""


 ############@ Annotate lineage & cell types logFC
 RUN_NAME=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
 #step 2: add lineage logFC stats
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --chr 1-22" >> ${JOBFILE}
echo ""

# step3: mash lineage logFC stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME}`
 echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME}"
 echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME}" >> ${JOBFILE}
 echo ""


##############
# step2: add celltype logFC stats
RUN_NAME_ANNOT=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22" >> ${JOBFILE}
echo ""


# step3: celltype logFC stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype" >> ${JOBFILE}
echo ""


##############
# step2: add lineage expr stats
RUN_NAME_ANNOT=lineage_condition___CellPropLineage_SVs_220409
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --chr 1-22" >> ${JOBFILE}
echo ""

# step3: lineage expr stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr" >> ${JOBFILE}
echo ""







##############
# step2: add celltype expr stats
RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr_celltype --chr`
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr_celltype --chr 1-22"
echo "${JOBID_AGG_1CHR} submitted ${OPTIONS} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr_celltype --chr 1-22" >> ${JOBFILE}
echo ""
# step3: celltype expr stats
JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} --dependency=afterok:${JOBID_AGG_1CHR} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R ${OPTIONS} --run_name ${RUN_NAME} --name expr_celltype`
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --name expr_celltype"
echo "${JOBID_AGG_MASH} submitted ${OPTIONS} --run_name ${RUN_NAME} --name expr_celltype" >> ${JOBFILE}
echo ""


  # SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
  # EQTL_SCRIPT_DIR="09_eQTLmapping"
  # EQTL_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/3_eQTL_mapping"
  #
  # USE_JOINT=FALSE
  # USE_CELLTYPE=TRUE




########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################



     ########@ run non conditional adjusted for cell type

########@ run non conditional adjusted for cell type
#
# RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}" >> ${JOBFILE}
# # echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
# RUN_NAME=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
# JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}`
# echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}" >> ${JOBFILE}
# echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype FAL"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
# RUN_NAME=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
# JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}`
# echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}" >> ${JOBFILE}
# echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
# # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE`
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE"
# # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE" >> ${JOBFILE}
# # echo ""
#
#
# ########
#
# RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
# SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
# # echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#  # JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --chr`
#  # echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT} "
#  # echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#  # echo ""
#
#  RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
#  JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr`
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT}"
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#  echo ""
#
#  RUN_NAME_ANNOT=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
#  JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr`
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} "
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} " >> ${JOBFILE}
#  echo ""
#
#  RUN_NAME_ANNOT=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
#  JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC_celltype --chr`
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} "
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} " >> ${JOBFILE}
#  echo ""
#
#
#  RUN_NAME=lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
#  SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
#  # echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}"
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#   echo ""
#
#   RUN_NAME_ANNOT=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --add_joint_mapping ${USE_JOINT} --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} "
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} " >> ${JOBFILE}
#   echo ""
#
#   RUN_NAME_ANNOT=lineage_condition___CellPropLineage_SVs_220409
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --add_joint_mapping ${USE_JOINT} --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --chr 1-22 --add_joint_mapping ${USE_JOINT}"
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#   echo ""
#
#   RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name expr_celltype --add_joint_mapping ${USE_JOINT} --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} "
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} " >> ${JOBFILE}
#   echo ""
#
#   RUN_NAME=celltype_condition___CellPropLineage_SVs_220409
#   SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
#   # echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#    JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --chr`
#    echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}"
#    echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#    echo ""
#
#
#   RUN_NAME=celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
#    SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
#    # echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#     JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --chr`
#     echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}"
#     echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT}" >> ${JOBFILE}
#     echo ""
#
#
#
# ######################
# # SBATCH_MASH_ARGS2=" --qos=geh --partition=gehbigmem -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_mash_FineMap_%A.log -J FineMap_mash"
# # --dependency=afterok:${JOBID_FINEMAP},${JOBID_MASH}"
# # echo "sbatch --mem=100G ${SBATCH_MASH_ARGS2}  ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c5_mashr_fineMapped.R --dist ${DIST} --run_name ${RUN_NAME}"
# # JOBID_MASH2=`sbatch --mem=100G ${SBATCH_MASH_ARGS2} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c5_mashr_fineMapped.R --dist ${DIST} --run_name ${RUN_NAME}`
# # echo "${JOBID_MASH2} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# # echo ""

# USE_CELLTYPE=TRUE
# USE_JOINT=FALSE
# RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
# SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR"
# # echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#  JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE} --chr`
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}" >> ${JOBFILE}
#  echo ""
#
#   RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --run_name_annotate ${RUN_NAME_ANNOT} --add_celltype_mapping ${USE_CELLTYPE} --name celltype --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}" >> ${JOBFILE}
#   echo ""
#
#   RUN_NAME_ANNOT=lineage_condition_logFC___CellPropLineage_SVs_220409
#   JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --run_name_annotate ${RUN_NAME_ANNOT} --add_celltype_mapping ${USE_CELLTYPE} --name logFC --chr`
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE} "
#   echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE} " >> ${JOBFILE}
#   echo ""
#

#
#
# RUN_NAME=lineage_condition___CellPropLineage_SVs_220409
# SBATCH_AGG_MASH_ARGS="--parsable --mem=100G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_mash_%A_%a.log -J eQTL_AGG_mash --dependency=afterok:$JOBID_AGG_1CHR"
# # echo "sbatch --array=1-22 ${SBATCH_AGG_MASH_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d6_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
#  JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}`
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}" >> ${JOBFILE}
#  echo ""
#
#  RUN_NAME_ANNOT=celltype_condition___CellPropLineage_SVs_220409
#  JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype`
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name celltype --chr 1-22 --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}" >> ${JOBFILE}
#  echo ""
#
#  RUN_NAME_ANNOT=lineage_condition_logFC___CellPropLineage_SVs_220409
#  JOBID_AGG_MASH=`sbatch ${SBATCH_AGG_MASH_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d7_mash_allCHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC`
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}"
#  echo "${JOBID_AGG_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME} --run_name_annotate ${RUN_NAME_ANNOT} --name logFC --add_joint_mapping ${USE_JOINT} --add_celltype_mapping ${USE_CELLTYPE}" >> ${JOBFILE}
#  echo ""
