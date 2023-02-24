#!/bin/bash
################################################################################
################################################################################
# File name: 3a8__launcher_part_a2_a3_a4.sh
# Author: Y.A., M.R., M.ON.
################################################################################
###############################################################################
# Step: run SuSiE and all chromosomes/celltypes/conditions and aggregate results
# (steps 3a2 to 3a4)
# same script is used for reQTLs
# Batch launcher script
################################################################################
################################################################################

################################################################################
#### Example calls
# ./3a8__launcher_part_a2_a3_a4.sh --dist 100000 --run_name lineage_condition___CellPropLineage_SVs_RUN1
# ./3a8__launcher_part_a2_a3_a4.sh --dist 100000 --run_name lineage_condition_logFC__CellPropLineage_SVs_RUN1
#
# ./3a8__launcher_part_a2_a3_a4.sh --dist 100000 --run_name celltype_condition___CellPropLineage_SVs_RUN1
# ./3a8__launcher_part_a2_a3_a4.sh --dist 100000 --run_name celltype_condition_logFC__CellPropLineage_SVs_RUN1

LOGS="./LOGS"
JOBFILE='./LOGS/launched_jobs_3a8.txt'

while [[ $# -gt 0 ]]; do
  case $1 in
    -r|--run_name)
      RUN_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    -d|--dist)
      DIST="$2"
      shift # past argument
      shift # past value
      ;;
    -j|--jobfile)
      JOBFILE="$2"
      shift # past argument
      shift # past value
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
    esac
  done

MISC_DIR='./MISC/'
EQTL_SCRIPT_DIR="./3__eQTL_mapping/3a_eQTL_mapping"
EQTL_DIR="3__eQTL_mapping/sumStats"

CELLSTATES=`ls ${EQTL_DIR}/${RUN_NAME} | grep '__'`
ALL_JOBID=''
#### run finemapping per chromosome & cellstate
for STATE in $CELLSTATES;
  do
SBATCH_FINEMAP_ARGS="--parsable --mem=150G --time=6:00:00 -o ${LOGS}/logfile_eQTL_FineMap_CTspecific_%A_%a.log -J eQTL_FineMap"
 echo "sbatch --array=1-22 ${SBATCH_FINEMAP_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a2_FineMapping_1CHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr \${SLURM_ARRAY_TASK_ID}" >> ${JOBFILE}
 JOBID_FINEMAP=`sbatch --array=1-22 ${SBATCH_FINEMAP_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a2_FineMapping_1CHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr`
 echo "${JOBID_FINEMAP} submitted --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr 1-22 "  >> ${JOBFILE}
 echo ""
  # compute FDR on each state
 SBATCH_FDR_ARGS="--parsable --mem=50G --time=4:00:00 -o ${LOGS}/logfile_eQTL_FDR_CTspecific_%A_%a.log -J eQTL_FDR --dependency=afterok:${JOBID_FINEMAP}"
 echo "sbatch --array=1-22 ${SBATCH_FDR_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a3_aggregate_FineMapping_allCHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE}">> ${JOBFILE}
 JOBID_FDR=`sbatch ${SBATCH_FDR_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/3a3__aggregate_FineMapping_allCHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE}`
 echo "${JOBID_FDR} submitted --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE}" >> ${JOBFILE}
 echo ""
 ALL_JOBID=${JOBID_FDR},${ALL_JOBID}
done;

# compile Finemapped eQTLs, extract 95% intervals and aggregate to identify independent eQTLs at each level
SBATCH_AGG_ARGS="--parsable --mem=60G -o ${LOGS}/logfile_eQTL_AGG_allCT_%A.log -J eQTL_AGG --dependency=afterok:${ALL_JOBID}"
 echo "sbatch ${SBATCH_AGG_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}"  >> ${JOBFILE}
 JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${MISC_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}`
 echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}" >> ${JOBFILE}
 echo ""

#####################################################################################################################
### after running 3a8_launcher for expression and logFC, at lineage & celltype level , go to 3a9_launcher.     ######
#####################################################################################################################
