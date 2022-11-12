#!/bin/bash
# not run
# sbatch --array=1-2222 --mem=80G --partition=common,dedicated -J MateQTL -o "log/logfile_MateQTLV2_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh

# sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage_SVs --run_id lineage_test --test TRUE
# sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage_SVs --logfc TRUE --run_id lineage_test --test TRUE
# sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage_SVs --run_id celltype_test --test TRUE

# 5 lineages logCPM
# sbatch --array=1-44 --mem=80G --partition=common --qos=normal -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage_SVs --run_id 220409
# 5 lineages logFC
# sbatch --array=1-44 --mem=80G --partition=common --qos=normal -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --state condition --covname lineage_condition__CellPropLineage_SVs --run_id 220409 --logfc TRUE
# 5 lineages logFC with logFC SVs
# sbatch --array=1-44 --mem=80G --partition=common --qos=normal -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --state condition --covname lineage_condition_logFC__CellPropLineage_SVs --run_id 220409 --logfc TRUE

# celltypes logCPM
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage_SVs --run_id 220409
# 22 celltypes, logFC
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage_SVs --run_id 220409 --logfc TRUE
# 22 celltypes, logFC with logFC SVs
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition_logFC__CellPropLineage_SVs --run_id 220409 --logfc TRUE

# celltypes logCPM pDC,cDC,ILC,plasma,...
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage_SVs --run_id 220409 --force_rerun FALSE
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --state condition --covname celltype_condition__CellPropLineage_SVs --run_id 220409 --logfc TRUE --force_rerun FALSE


# ## test RUN, only lauch script is run
# sbatch --array=1-44 --mem=4G --time=00:59:00 --partition=common,geh --qos=fast -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --cellprop TRUE --test TRUE

TASK_NB=${SLURM_ARRAY_TASK_ID}
# TASK_NB=1
# for TASK_NB in `seq 1 2222`;
# do
PERM=`echo \($TASK_NB-1\)/22 | bc`
CHR=`echo $TASK_NB-22*\($PERM\) | bc`
echo "TASK_NB=$TASK_NB, PERM=$PERM, CHR=$CHR"
# done

SCRATCH='/pasteur/appa/scratch/mrotival/sc_Process/'
mkdir $SCRATCH
WORK_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/"
OUT_DIR="$WORK_DIR/project/pop_eQTL/data/3_eQTL_mapping"
SCRIPT_DIR="${WORK_DIR}/resources/template_scripts/processing_pipeline/09_eQTLmapping"

CELLTYPE="celltype"
STATE='condition'
COV_RUN_NAME='celltype_condition__CellPropLineage_SVs'
RUN_ID='test'
LOGFC='FALSE'
FORCE='TRUE' # should results be recomputed if already present ?

### read command line arguments to define $CELLTYPE, $STATE, $RUN_ID, and $COV_RUN_NAME
echo "reading args"
CMD=""

while [[ $# -gt 0 ]]; do
  echo "'\$1 \$2' is $1 $2"
  case $1 in
    -t|--celltype)
    CELLTYPE="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    -a|--state)
    STATE="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    -r|--run_id)
    RUN_ID="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    -v|--covname)
    COV_RUN_NAME="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    -w|--workdir)
    OUT_DIR="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    -f|--logfc)
    LOGFC="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
    --force_rerun)
    FORCE="$2"
    CMD="${CMD} ${1} ${2}"
    shift
    shift
    ;;
     *)
    CMD="${CMD} ${1}"
    shift
  esac
done

echo "$CMD"

### define RUN_NAME
RUN_NAME="${COV_RUN_NAME}_${RUN_ID}"
COV_RUN_NAME_SIMPLIFIED=$(echo ${COV_RUN_NAME} | sed -e "s/"${CELLTYPE}"\_"${STATE}"//")

echo ${RUN_NAME}

# define condition columns for the cell metadata file
# (condition is named COND in the cell metadata file)
if [ ${STATE} = "condition" ]; then
  STATE_COL='COND'
else
  STATE_COL=${STATE}
fi

META_CELLS="${WORK_DIR}/project/pop_eQTL/data/2_population_differences/sce_clean_metadata.tsv"
echo ${META_CELLS} ${STATE_COL} ${CELLTYPE}
### define the set of celltypes and states to Loop on based on ${CELLTYPE} & ${STATE}
#(we remove the colnames and mDC lines)
echo 'listing cellstates to test'
#awk -V

awk -F'\t' -v STATE=${STATE_COL} -v CELLTYPE=${CELLTYPE} 'BEGIN { OFS = "\t" }
                                    (NR==1){for(i=1;i<=NF;i++){
                                      if($i==CELLTYPE) celltype=i;
                                      if($i==STATE) state=i;
                                      }
                                    }{print $(celltype),$(state)}' $META_CELLS | grep -v 'celltype' | grep -v 'OTHER' | grep -v 'lineage' | sort | uniq > ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_prov.txt

# remove cell types that have no covariates
# cat ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_prov.txt | grep -v 'ILC' | grep -v 'pDC' | grep -v 'cDC' | grep -v 'Plasmablast' | grep -v 'Progenitor'> ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

cat ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_prov.txt > ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

echo 'listing complete'
cat ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
echo 'printed'

# load required modules
module load R/4.1.0

while read line; do
  myCELLTYPE=`echo "$line" | cut -f 1`
  mySTATE=`echo "$line" | cut -f 2`
  echo "$myCELLTYPE - $mySTATE"
  if [[ "$LOGFC" == "TRUE" ]]; then LOGFC_FLAG='_logFC'; else LOGFC_FLAG=''; fi
  TARGET_FILE="${OUT_DIR}/${CELLTYPE}_${STATE}${LOGFC_FLAG}___${COV_RUN_NAME_SIMPLIFIED}_${RUN_ID}/${myCELLTYPE}__${mySTATE}/perm/eQTL_ALL_chr${CHR}_bestP_perFeature_perm${PERM}.txt.gz"
  echo "searching for $TARGET_FILE"
  if [[ "$LOGFC" == "TRUE" && "$mySTATE" == "NS" ]]; then echo "skipping..."; fi
  if [[ "$LOGFC" == "TRUE" && "$mySTATE" == "NS" ]]; then continue; fi
  if [[ -f $TARGET_FILE && "$FORCE" != "TRUE" ]]; then continue; fi
  echo "Rscript ${SCRIPT_DIR}/09b_MatrixEQTL_core_V4.R --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --perm ${PERM} --chr ${CHR} ${CMD}" >> ${SCRATCH}/myJobs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
#  echo "Rscript ${SCRIPT_DIR}/09b_MatrixEQTL_core_V4.R --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --perm ${PERM} --chr ${CHR} ${CMD}"
  Rscript ${SCRIPT_DIR}/09b_MatrixEQTL_core_V4.R --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --perm ${PERM} --chr ${CHR} ${CMD}
  # Rscript ${SCRIPT_DIR}/09b_MatrixEQTL_core_V4.R --celltype ${CELLTYPE} --state ${STATE} --covdir ${COVDIR} --covname ${COVNAME} --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --chr ${CHR} --perm ${PERM} --dist ${DIST} --run_id ${RUN_ID} --cellprop ${ADD_CELLPROP} --test ${TEST} --workdir ${WORKDIR} --scratch ${SCRATCH} --expression ${EXPRESSION_FILE} --metadata ${METADATA_FILE}

done < ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

chmod -R 775 ${OUT_DIR}/${CELLTYPE}_${STATE}${LOGFC_FLAG}___${COV_RUN_NAME_SIMPLIFIED}_${RUN_ID}

#chmod -R 775 ${OUT_DIR}/${RUN_NAME}
# rm ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

exit



# CLUST=0
# CHR=1
# PERM=0
# time Rscript ${SCRIPT_DIR}/MatrixEQTL_core.R --clust $CLUST --chr ${CHR} --perm ${PERM} --dist 100000



############@ submitted jobs 090222_test
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --run_id 090222_test --test TRUE
# Submitted batch job 64865067
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_noSV --cellprop TRUE --run_id 090222_test --test TRUE
# Submitted batch job 64865068
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --cellprop TRUE --run_id 090222_test --test TRUE
# Submitted batch job 64865069
# [mrotival@maestro-submit processing_pipeline]$
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --run_id 090222_logfc_test --test TRUE --logfc TRUE
# Submitted batch job 64865070
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222_test --test TRUE
# Submitted batch job 64865071
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1,22,23,44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222_test --test TRUE --logfc TRUE
# Submitted batch job 64865073



############@ submitted jobs 090222 actual run
# [mrotival@maestro-submit processing_pipeline]$ ###### 8 celltypes logFC
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --run_id 090222 --logfc TRUE
# Submitted batch job 64893627
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_noSV --cellprop TRUE --run_id 090222 --logfc TRUE
# Submitted batch job 64893628
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --cellprop TRUE --run_id 090222 --logfc TRUE
# Submitted batch job 64893629
# [mrotival@maestro-submit processing_pipeline]$
# [mrotival@maestro-submit processing_pipeline]$ ###### 19 celltypes log CPM
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222
# Submitted batch job 64893630
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_noSV --cellprop TRUE --run_id 090222
# Submitted batch job 64893631
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --cellprop TRUE --run_id 090222
# Submitted batch job 64893632
# [mrotival@maestro-submit processing_pipeline]$
# [mrotival@maestro-submit processing_pipeline]$ ###### 19 celltypes log FC
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222 --logfc TRUE
# Submitted batch job 64893633
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_noSV --cellprop TRUE --run_id 090222 --logfc TRUE
# Submitted batch job 64893634
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --cellprop TRUE --run_id 090222 --logfc TRUE
# Submitted batch job 64893635


#### rerun, no cell prop adjustment
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --run_id 090222
# Submitted batch job 65406385
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.8level --state condition --covname Covariate_CT8_condition_nosubset_noprop_proptype_defSV --run_id 090222 --logfc TRUE
# Submitted batch job 65406377
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222
# Submitted batch job 65406378
# [mrotival@maestro-submit processing_pipeline]$ sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV3_%A_%a.log" ./09a_run_eQTL_MatrixEQTL_V3.sh --celltype celltype.19level --state condition --covname Covariate_CT19_condition_nosubset_noprop_proptype_defSV --run_id 090222 --logfc TRUE
# Submitted batch job 65406379
