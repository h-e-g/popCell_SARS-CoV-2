#!/bin/bash
# ./09d4_launcher.sh --dist 100000 --run_name lineage_condition___CellPropLineage_SVs_220409
# ./09d4_launcher.sh --dist 100000 --run_name lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
#
# ./09d4_launcher.sh --dist 100000 --run_name celltype_condition___CellPropLineage_SVs_220409
# ./09d4_launcher.sh --dist 100000 --run_name celltype_condition_logFC__logFC__CellPropLineage_SVs_220409


# retired (adjusted on wrong SVs)
##old# ./09d4_launcher.sh --dist 100000 --run_name lineage_condition_logFC___CellPropLineage_SVs_220409
##old# ./09d4_launcher.sh --dist 100000 --run_name celltype_condition_logFC___CellPropLineage_SVs_220409

# previous run 240322
JOBFILE='./zz_launchedjobs_20052022.txt'

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

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
EQTL_SCRIPT_DIR="09_eQTLmapping"
EQTL_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/3_eQTL_mapping"

CELLSTATES=`ls $EQTL_DIR/$RUN_NAME | grep '__'`
ALL_JOBID=''
#### run finemapping per chromosome & cellstate
for STATE in $CELLSTATES;
  do
SBATCH_FINEMAP_ARGS="--parsable --mem=150G --qos=normal --partition=common --time=6:00:00 -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_FineMap_CTspecific_%A_%a.log -J eQTL_FineMap"
 echo "sbatch --array=1-22 ${SBATCH_FINEMAP_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d1_FineMapping_1CHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr \${SLURM_ARRAY_TASK_ID}" >> ${JOBFILE}
 JOBID_FINEMAP=`sbatch --array=1-22 ${SBATCH_FINEMAP_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d1_FineMapping_1CHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr`
 echo "${JOBID_FINEMAP} submitted --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --chr 1-22 "  >> ${JOBFILE}
 echo ""
  # compute FDR on each state
 SBATCH_FDR_ARGS="--parsable --mem=50G --qos=geh --partition=geh --time=4:00:00 -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_FDR_CTspecific_%A.log -J eQTL_FDR --dependency=afterok:${JOBID_FINEMAP}"
 echo "sbatch --array=1-22 ${SBATCH_FDR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d2_aggregate_allCHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE} --dependency=afterok:${JOBID_FINEMAP}">> ${JOBFILE}
 JOBID_FDR=`sbatch ${SBATCH_FDR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d2_aggregate_allCHR_1cond.R --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE}`
 echo "${JOBID_FDR} submitted --dist ${DIST} --run_name ${RUN_NAME} --cellstate ${STATE}" >> ${JOBFILE}
 echo ""
# ALL_JOBID=${JOBID_FDR},${ALL_JOBID}
done;

# compile Finemapped eQTLs, extract 95% intervals and aggregate to identify independent eQTLs at each level
SBATCH_AGG_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_%A.log -J eQTL_AGG --dependency=afterok:${JOBID_FDR}"
 echo "sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}"  >> ${JOBFILE}
 JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME}`
 echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}"
 echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME}" >> ${JOBFILE}
 echo ""

#####################################################################################################################
### after running 09d4_launcher for expression and logFC, at lineage & celltype level , go to 09d8lancher.     ######
#####################################################################################################################






 # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --correct_celltype TRUE`
 #  echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE"
 # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional FALSE --correct_celltype TRUE" >> ${JOBFILE}
 # # echo ""
 # JOBID_AGG=`sbatch ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d3_aggregate_allCHR_allcond.R --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE`
 # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE"
 # echo "${JOBID_AGG} submitted --dist ${DIST} --run_name ${RUN_NAME} --conditional TRUE --correct_celltype TRUE" >> ${JOBFILE}
 # echo ""

# SBATCH_AGG_1CHR_ARGS="--parsable --mem=60G --qos=geh --partition=geh -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_AGG_allCT_1Chr_%A_%a.log -J eQTL_AGG_1CHR --dependency=afterok:${JOBID_AGG}"
#  echo "sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d4_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}" >> ${JOBFILE}
#  JOBID_AGG_1CHR=`sbatch --array=1-22 ${SBATCH_AGG_1CHR_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09d4_aggregate_1CHR_allCond.R --dist ${DIST} --run_name ${RUN_NAME} --chr`
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 "
#  echo "${JOBID_AGG_1CHR} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22 " >> ${JOBFILE}
#  echo ""

#### run mash on eQTL best SNP
 # SBATCH_MASH_ARGS="--parsable --mem=200G --qos=geh --partition=gehbigmem -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_mash_%A.log -J eQTL_mash --dependency=afterok:${JOBID}"
 # echo "sbatch ${SBATCH_MASH_ARGS}  ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c3_mash_allCHR.R --dist ${DIST} --run_name ${RUN_NAME}"
 # JOBID_MASH=`sbatch ${SBATCH_MASH_ARGS}  ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c3_mash_allCHR.R --dist ${DIST} --run_name ${RUN_NAME}`
 # echo "${JOBID_MASH} submitted --dist ${DIST} --run_name ${RUN_NAME}"
 # echo ""

#### run mash on finemapped eQTLs
# SBATCH_MASH_ARGS2=" --qos=geh --partition=gehbigmem -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_mash_FineMap_%A.log -J FineMap_mash"
# --dependency=afterok:${JOBID_FINEMAP},${JOBID_MASH}"
# echo "sbatch --mem=100G ${SBATCH_MASH_ARGS2}  ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c5_mashr_fineMapped.R --dist ${DIST} --run_name ${RUN_NAME}"
# JOBID_MASH2=`sbatch --mem=100G ${SBATCH_MASH_ARGS2} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c5_mashr_fineMapped.R --dist ${DIST} --run_name ${RUN_NAME}`
# echo "${JOBID_MASH2} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# echo ""

 #### annotate eQTL best SNPs with reQTL stats
# SBATCH_REQTL_ARGS="--parsable --qos=geh --partition=gehbigmem -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_mash_reQTL_%A.log -J reQTL_mash --dependency=afterok:${JOBID_MASH}"
# echo "sbatch --mem=100G ${SBATCH_REQTL_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09e_add_logFC.R --dist ${DIST} --run_name ${RUN_NAME}"
# JOBID_REQTL=`sbatch --mem=100G ${SBATCH_REQTL_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09e_add_logFC.R --dist ${DIST} --run_name ${RUN_NAME}`
# echo "${JOBID_REQTL} submitted --dist ${DIST} --run_name ${RUN_NAME}"
# echo ""

# #### annotate finemapped eQTL with reQTL stats
# SBATCH_REQTL_ARGS="--parsable --qos=geh --partition=gehbigmem -o ${EQTL_SCRIPT_DIR}/log/logfile_eQTL_mash_reQTL_%A.log -J reQTL_mash --dependency=afterok:${JOBID_MASH}"
# echo "sbatch --mem=100G ${SBATCH_REQTL_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09e_add_logFC.R --dist ${DIST} --run_name ${RUN_NAME}"
# JOBID_REQTL=`sbatch --mem=100G ${SBATCH_REQTL_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09e_add_logFC.R --dist ${DIST} --run_name ${RUN_NAME}``
# echo "${JOBID_REQTL} submitted --dist ${DIST} --run_name ${RUN_NAME}"


# sbatch --array=1-22 --qos=ultrafast --partition=common,dedicated,geh,gehbigmem -o log/test.log -J test ./00_Rscript.sh Z_printargs.R --dist 10000 --run_name 'test_run' --chr
