#!/bin/bash
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END
#
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --run_name lineage_condition___CellPropLineage_SVs_220409 --dist 100000
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name celltype_condition___CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name celltype_condition_logFC__logFC__CellPropLineage_SVs_220409

# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name lineage_condition___CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name celltype_condition___CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 100000 --run_name celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
# # job IDS= 27832975-27832977, 27833085
#
# SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
# EQTL_SCRIPT_DIR="09_eQTLmapping"
# sbatch --array=2-946 --parsable --mem=80G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 250000 --run_name lineage_condition___CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=80G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 250000 --run_name celltype_condition___CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=80G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 250000 --run_name celltype_condition_logFC__logFC__CellPropLineage_SVs_220409
# sbatch --array=2-946 --parsable --mem=80G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_coloc_%A_%a.log -J coloc ./11c3_COVIDloci_coloc_launcher.sh --dist 250000 --run_name lineage_condition_logFC__logFC__CellPropLineage_SVs_220409
# job IDS= 27833481-27833484

#TODO: look for Errors (TIME LIMIT | memory) in logs

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

SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
EQTL_SCRIPT_DIR="09_eQTLmapping"
SCRATCH='/pasteur/appa/scratch/mrotival/sc_Process/'

# load required modules
module load R/4.1.0
module load samtools
EQTL_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/3_eQTL_mapping"
FILE="${EQTL_DIR}/overlap_COVID_1e5_within_100kb_of_eQTL_and_reQTL_peak_snp_gene_only.txt"

LINE=`head -n ${SLURM_ARRAY_TASK_ID} $FILE | tail -n 1`
#while read line; do
  #myCELLTYPE=`echo "$line" | cut -f 1`
  #mySTATE=`echo "$line" | cut -f 2`
  mySNP=`echo "$LINE" | cut -f 1`
  myGENE=`echo "$LINE" | cut -f 2`
  mySymbol=`echo "$LINE" | cut -f 3`
  myCHR=`echo "$LINE" | cut -f 4`

  echo "$mySymbol - $mySNP - $myCHR"
  echo "Rscript ${SCRIPT_DIR}/11c2_COVIDloci_coloc_1eQTL.R --snp ${mySNP} --gene ${myGENE} --name ${mySymbol} --chr ${myCHR} --run_name ${RUN_NAME}" >> ${SCRATCH}/myJobs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
  Rscript ${SCRIPT_DIR}/11c2_COVIDloci_coloc_1eQTL.R --snp ${mySNP} --gene ${myGENE} --name ${mySymbol} --chr ${myCHR} --run_name ${RUN_NAME} --dist $DIST


  # R ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/11c2_COVIDloci_coloc_1eQTL.R --snp rs114839171 --gene ENSG00000246451 --name AL049840.2 --chr 14 --run_name celltype_condition___CellPropLineage_SVs_220409 --dist 100000
#done

#### run aggregate per chromosome
# SBATCH_AGG_ARGS="--parsable --mem=60G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_eQTL_agg_%A_%a.log -J eQTL_ag${CHR}"
# echo "sbatch --array=1-22 ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c1_MatrixEQTL_aggregate_1CHR.R --dist ${DIST} --run_name ${RUN_NAME} --chr \${SLURM_ARRAY_TASK_ID}"
# JOBID=`sbatch --array=1-22 ${SBATCH_AGG_ARGS} ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/09c1_MatrixEQTL_aggregate_1CHR.R --dist ${DIST} --run_name ${RUN_NAME} --chr`
# echo "${JOBID} submitted --dist ${DIST} --run_name ${RUN_NAME} --chr 1-22"
# echo ""
