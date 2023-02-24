#!/bin/bash
################################################################################
################################################################################
# File name: 6b_COVIDloci_coloc_launcher.sh
# Author:  Y.A., M.R., M.ON., J.MR.
################################################################################
################################################################################
# Step: extract sumstats and test colocalization for one eQTL
# launcher script
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=6b
#SBATCH --output=${LOG_DIR}/6b__%J.log
#SBATCH -J coloc_1eQTL
#SBATCH --mem=60G
################################################################################
################################################################################
# Command history
# sbatch --array=2-946 ./6b_COVIDloci_coloc_launcher.sh --dist 250000 --run_name lineage_condition___CellPropLineage_SVs_RUN1
# sbatch --array=2-946 ./6b_COVIDloci_coloc_launcher.sh --dist 250000 --run_name lineage_condition_logFC__CellPropLineage_SVs_RUN1
#
# sbatch --array=2-946 ./6b_COVIDloci_coloc_launcher.sh --dist 250000 --run_name celltype_condition___CellPropLineage_SVs_RUN1
# sbatch --array=2-946 ./6b_COVIDloci_coloc_launcher.sh --dist 250000 --run_name celltype_condition_logFC__CellPropLineage_SVs_RUN1

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
   *)
     POSITIONAL_ARGS+=("$1") # save positional arg
     shift # past argument
     ;;
   esac
 done

# load required modules
module load R/4.1.0
module load samtools
FILE="6__COVID19_association/data/overlap_reQTL_and_reQTL__COVID_loci_P1e5__100kb.txt"

LINE=`head -n ${SLURM_ARRAY_TASK_ID} $FILE | tail -n 1`
while read line; do
  mySNP=`echo "$LINE" | cut -f 1`
  myGENE=`echo "$LINE" | cut -f 2`
  mySymbol=`echo "$LINE" | cut -f 3`
  myCHR=`echo "$LINE" | cut -f 4`

  echo "$mySymbol - $mySNP - $myCHR"
  echo "Rscript ./6b_COVIDloci_coloc_1eQTL.R --snp ${mySNP} --gene ${myGENE} --name ${mySymbol} --chr ${myCHR} --run_name ${RUN_NAME}"
  Rscript ${SCRIPT_DIR}/6b_COVIDloci_coloc_1eQTL.R --snp ${mySNP} --gene ${myGENE} --name ${mySymbol} --chr ${myCHR} --run_name ${RUN_NAME} --dist $DIST

done
