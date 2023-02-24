#!/bin/bash
################################################################################
################################################################################
# File name: 3a1__get_eQTL_SumStats_MatrixEQTL.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use MatrixEQTL models to estimate eQTL effect sizes, R2, and p-value
# same script is used for reQTLs
# Batch launcher script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
SCRATCH='./SCRATCH/'
mkdir $SCRATCH
################################################################################
# SLURM options
#SBATCH --job-name=3a1
#SBATCH --output=${LOG_DIR}/3a1__%J.log
#SBATCH --mem=80G
#SBATCH --array=1-44
################################################################################

################################################################################
# Setup
# load modules
module load R/4.1.0

################################################################################
# Command history
# 5 lineages logCPM
# sbatch ./3a1__get_eQTL_SumStats_MatrixEQTL.sh --celltype lineage --covname lineage_condition___CellPropLineage_SVs --run_id RUN1
# 5 lineages logFC
# sbatch ./3a1__get_eQTL_SumStats_MatrixEQTL.sh --celltype lineage --covname lineage_condition_logFC__CellPropLineage_SVs --run_id RUN1 --logfc TRUE

# celltypes logCPM
# sbatch ./3a1__get_eQTL_SumStats_MatrixEQTL.sh --celltype celltype --covname celltype_condition___CellPropLineage_SVs --run_id RUN1
# 22 celltypes, logFC
# sbatch ./3a1__get_eQTL_SumStats_MatrixEQTL.sh --celltype celltype --covname celltype_condition_logFC__CellPropLineage_SVs --run_id RUN1 --logfc TRUE

# TASK NB corresponds to CHR+ 22*PERM
TASK_NB=${SLURM_ARRAY_TASK_ID}
# obtain PERMutation number and CHRomosome
PERM=`echo \($TASK_NB-1\)/22 | bc`
CHR=`echo $TASK_NB-22*\($PERM\) | bc`
echo "TASK_NB=$TASK_NB, PERM=$PERM, CHR=$CHR"

OUT_DIR="3__eQTL_mapping"

###fixed parameter values:
# column featuring the sample condition in sce_clean_metadata
STATE_COL='COND'

### arg parameter values:
# set default parameter values
CELLTYPE="celltype" # should eQTL be run at celltype (22 celltype) or lineage level (5 lineages)
COV_RUN_NAME='lineage_condition__CellPropLineage_SVs' # folder where covariates can be found
RUN_ID='RUN1' # ID of the run (to avoid confusion)
LOGFC='FALSE' # Should eQTL be run on log fold change to NS ? (TRUE or FALSE)

### read command line arguments to define $CELLTYPE, $RUN_ID, $LOGFC and $COV_RUN_NAME
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
echo ${RUN_NAME}

# define condition columns for the cell metadata file
# (condition is named COND in the cell metadata file)

#
##### use cell meta data defined in 1b1__celltype_identification.R to deifne the set of cell types to consider
META_CELLS="${AGGR_QC_DIR}/data/sce_clean_metadata.tsv"
### define the set of celltypes and states to Loop on based on ${CELLTYPE} & ${STATE}
echo 'listing cellstates to test'


awk -F'\t' -v STATE=${STATE_COL} -v CELLTYPE=${CELLTYPE} 'BEGIN { OFS = "\t" }
                                    (NR==1){for(i=1;i<=NF;i++){
                                      if($i==CELLTYPE) celltype=i;
                                      if($i==STATE) state=i;
                                      }
                                    }{print $(celltype),$(state)}' $META_CELLS | grep -v 'celltype' | grep -v 'lineage' | sort | uniq > ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

echo 'listing complete'
cat ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
echo 'printed'

# load required modules
module load R/4.1.0

while read line; do
  myCELLTYPE=`echo "$line" | cut -f 1`
  mySTATE=`echo "$line" | cut -f 2`
  echo "$myCELLTYPE - $mySTATE"

  Rscript ./3a1__get_eQTL_SumStats_MatrixEQTL.R --mycelltype ${myCELLTYPE} --mystate ${mySTATE} --perm ${PERM} --chr ${CHR} ${CMD}

done < ${SCRATCH}/CLUSTER_LIST_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

rm -r $SCRATCH

exit
