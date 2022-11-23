#!/bin/bash
################################################################################
################################################################################
# File name: 0a1__STARsolo_align.sh
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Per-library alignment of FASTQ reads, valid barcode detection
# Effector script
################################################################################
################################################################################

LOG_DIR="../../../LOG"
LIB_COUNT=135

################################################################################
# SLURM options
#SBATCH --job-name=0a1
#SBATCH --output=${LOG_DIR}/0a1__%A_%a.log
#SBATCH --mem=140000
#SBATCH --ntasks=16
#SBATCH --array=1-${LIB_COUNT}
################################################################################

LIB=${SLURM_ARRAY_TASK_ID}

################################################################################
# Setup

# load modules
module load STAR/2.7.9a

PROJECT="HN00163853"
EVO_IMMUNO_POP="/pasteur/zeus/projets/p02/evo_immuno_pop"
RAWDIR="${EVO_IMMUNO_POP}/single_cell/raw"
FASTQDIR="${RAWDIR}/${PROJECT}/${PROJECT}_10X_RawData_Outs"
REFDIR="${EVO_IMMUNO_POP}/single_cell/resources/references/RNA/update_2021/human_iav_sars/star/"
OUTDIR="${EVO_IMMUNO_POP}/single_cell/aligned"

# run once

#cd ${RAWDIR}
#ls -l | awk -F" " '{print $9}' | grep -e "^R[1-4]" > scrnaseq_sample_list.txt

# get lane and flowcell identifiers per sample (i.e., library)

SAMPLE=`head -n ${SAMPLE_NUM} ${FASTQDIR}/scrnaseq_sample_list.txt | tail -n 1 | cut -d ' ' -f 1`

SAMPLEDIR=${FASTQDIR}/${SAMPLE}

rm ${SAMPLEDIR}/${SAMPLE}_fcell_list.txt

ls -l ${SAMPLEDIR} | awk -F" " '{print $9}' | tail -n -2 | grep -v "_fcell_list\.txt" | grep -e "^." > ${SAMPLEDIR}/${SAMPLE}_fcell_list.txt

NFLOWCELL=`wc -l ${SAMPLEDIR}/${SAMPLE}_fcell_list.txt | cut -d" " -f1`

if [ $NFLOWCELL -gt 1 ]
then

	declare -a ARRAY_FLOWCELL
	for i in $(seq 1 $NFLOWCELL)
	do
		FLOWCELL=`head -n ${i} ${SAMPLEDIR}/${SAMPLE}_fcell_list.txt | tail -n 1 | cut -d" " -f1`
		ARRAY_FLOWCELL[${i}]="${SAMPLEDIR}/${FLOWCELL}"
	done

	printf -v FLOWCELLPATH '%s,' "${ARRAY_FLOWCELL[@]}"
	FLOWCELLPATH="${FLOWCELLPATH::-1}"

elif [ $NFLOWCELL -eq 1 ]
then
	FLOWCELL=`head -n 1 ${SAMPLEDIR}/${SAMPLE}_fcell_list.txt`
	FLOWCELLPATH="${SAMPLEDIR}/${FLOWCELL}"
	SNO=`ls $SAMPLEDIR/$FLOWCELL | grep -e "${SAMPLE}" | head -n 1 | sed 's/_/\n/g' | grep -E "^S[0-9]"`
	LIBID=`echo $SAMPLE | cut -d"_" -f1`
fi

################################################################################
# Run STARsolo: alignment

echo "Project: ${PROJECT}"
echo "Sample: ${SAMPLE}"

start=`date +%s`

STAR --genomeDir ${REFDIR} \
	--readFilesCommand zcat \
	--readFilesIn ${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L001_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L002_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L003_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L004_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L005_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L006_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L007_R2_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L008_R2_001.fastq.gz ${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L001_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L002_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L003_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L004_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L005_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L006_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L007_R1_001.fastq.gz,${FASTQDIR}/${SAMPLE}/${FLOWCELL}/${SAMPLE}_${SNO}_L008_R1_001.fastq.gz \
	--outFileNamePrefix ${OUTDIR}/${PROJECT}/${LIBID}. \
	--runThreadN 16 \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist ${EVO_IMMUNO_POP}/single_cell/resources/important_files/3M-february-2018.txt \
	--soloUMIlen 12 \
	--clipAdapterType CellRanger4 \
	--outFilterScoreMin 30 \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	--soloUMIfiltering MultiGeneUMI_CR \
	--soloUMIdedup 1MM_CR \
	--soloCellFilter None \
	--soloFeatures Gene GeneFull \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	--outSAMtype BAM SortedByCoordinate

chmod -R 775 ${OUTDIR}/${PROJECT}/${LIBID}.Solo.out
chmod -R 775 ${OUTDIR}/${PROJECT}/${LIBID}.Aligned.sortedByCoord.out.bam

end=`date +%s`
runtime=$((end-start))
echo 'finished alignment'
echo $runtime

################################################################################
# Run STARsolo: cell filtering

# load previous version of STAR without bug in cell filtering
module purge
module load STAR/2.7.8a

start=`date +%s`

STAR --runMode soloCellFiltering ${OUTDIR}/${PROJECT}/${LIBID}.Solo.out/Gene/raw/ \
	${OUTDIR}/${PROJECT}/${LIBID}.Solo.out/Gene/filtered/ \
	--soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000

mv ${OUTDIR}/${PROJECT}/${LIBID}.Solo.out/Gene/filtered/features.tsv ${OUTDIR}/${PROJECT}/${LIBID}.Solo.out/Gene/filtered/genes.tsv
chmod -R 775 ${OUTDIR}/${PROJECT}/${LIBID}.Solo.out/Gene/filtered

end=`date +%s`
runtime=$((end-start))
echo 'finished empty drops cell identification'
echo $runtime
