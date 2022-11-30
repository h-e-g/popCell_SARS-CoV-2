#!/bin/bash

# 5 lineages logCPM
# sbatch --array=1-44 --mem=80G --partition=common --qos=normal -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./3a1_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --covname lineage_condition__CellPropLineage_SVs --run_id 220409
# 5 lineages logFC
# sbatch --array=1-44 --mem=80G --partition=common --qos=normal -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./3a1_run_eQTL_MatrixEQTL_V4.sh --celltype lineage --covname lineage_condition_logFC__CellPropLineage_SVs --run_id 220409 --logfc TRUE

# celltypes logCPM
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./3a1_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --covname celltype_condition__CellPropLineage_SVs --run_id 220409
# 22 celltypes, logFC
# sbatch --array=1-44 --mem=80G --partition=geh --qos=geh -J MateQTL -o "log/logfile_MateQTLV4_%A_%a.log" ./3a1_run_eQTL_MatrixEQTL_V4.sh --celltype celltype --covname celltype_condition_logFC__CellPropLineage_SVs --run_id 220409 --logfc TRUE

# TASK NB corresponds to CHR+ 22*PERM
TASK_NB=${SLURM_ARRAY_TASK_ID}

PERM=`echo \($TASK_NB-1\)/22 | bc`
CHR=`echo $TASK_NB-22*\($PERM\) | bc`
echo "TASK_NB=$TASK_NB, PERM=$PERM, CHR=$CHR"

SCRATCH='/pasteur/appa/scratch/mrotival/sc_Process/'
mkdir $SCRATCH
WORK_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/"
OUT_DIR="$WORK_DIR/project/pop_eQTL/data/3_eQTL_mapping"
SCRIPT_DIR="${WORK_DIR}/resources/template_scripts/processing_pipeline/09_eQTLmapping"


# set default parameter values
CELLTYPE="celltype" # should eQTL be run at celltype (22 celltype) or lineage level (5 lineages)
COV_RUN_NAME='celltype_condition__CellPropLineage_SVs' # folder where covariates can be found
RUN_ID='RUN1' # ID of the run (to avoid confusion)
LOGFC='FALSE' # Should eQTL be run on log fold change to NS ? (TRUE or FALSE)

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
     *)
    CMD="${CMD} ${1}"
    shift
  esac
done

echo "$CMD"

### define RUN_NAME
RUN_NAME="${COV_RUN_NAME}_${RUN_ID}"
COV_RUN_NAME_SIMPLIFIED=$(echo ${COV_RUN_NAME} | sed -e "s/"${CELLTYPE}"\_condition//")

echo ${RUN_NAME}

# define condition columns for the cell metadata file
# (condition is named COND in the cell metadata file)

META_CELLS="${WORK_DIR}/project/pop_eQTL/data/2_population_differences/sce_clean_metadata.tsv"
echo ${META_CELLS} ${STATE_COL} ${CELLTYPE}
### define the set of celltypes and states to Loop on based on ${CELLTYPE} & ${STATE}
#(we remove the colnames and mDC lines)
echo 'listing cellstates to test'
#awk -V

STATE_COL='COND'

awk -F'\t' -v STATE=${STATE_COL} -v CELLTYPE=${CELLTYPE} 'BEGIN { OFS = "\t" }
                                    (NR==1){for(i=1;i<=NF;i++){
                                      if($i==CELLTYPE) celltype=i;
                                      if($i==STATE) state=i;
                                      }
                                    }{print $(celltype),$(state)}' $META_CELLS | grep -v 'celltype' | grep -v 'OTHER' | grep -v 'lineage' | sort | uniq > /CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

echo 'listing complete'
cat ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
echo 'printed'

# load required modules
module load R/4.1.0

while read line; do
  myCELLTYPE=`echo "$line" | cut -f 1`
  mySTATE=`echo "$line" | cut -f 2`
  echo "$myCELLTYPE - $mySTATE"
  if [[ "$LOGFC" == "TRUE" ]];
    then LOGFC_FLAG='_logFC';
    else LOGFC_FLAG='';
  fi

  TARGET_FILE="${OUT_DIR}/${CELLTYPE}_${STATE}${LOGFC_FLAG}___${COV_RUN_NAME_SIMPLIFIED}_${RUN_ID}/${myCELLTYPE}__${mySTATE}/perm/
  eQTL_ALL_chr${CHR}_bestP_perFeature_perm${PERM}.txt.gz"

  echo "searching for $TARGET_FILE"
  if [[ "$LOGFC" == "TRUE" && "$mySTATE" == "NS" ]];
    then echo "skipping...";
         continue;
  fi

  if [[ -f $TARGET_FILE && "$FORCE" != "TRUE" ]];
    then continue;
  fi

  Rscript ${SCRIPT_DIR}/3a1_run_eQTL_MatrixEQTL.R --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --perm ${PERM} --chr ${CHR} ${CMD}

done < ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

chmod -R 775 ${OUT_DIR}/${CELLTYPE}_${STATE}${LOGFC_FLAG}___${COV_RUN_NAME_SIMPLIFIED}_${RUN_ID}
exit
