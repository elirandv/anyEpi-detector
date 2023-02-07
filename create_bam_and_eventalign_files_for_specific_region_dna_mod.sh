EXPERIMENT_ROOT_DIR=$1
REF_SEQ=$2
PKL_FILE_MODEL=$3
BARCODE=$4
CHR=$5
START=$6
END=$7

WORKINK_DIR_PREFIX="${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}"
FAST5_FILES_DIR="${WORKINK_DIR_PREFIX}/"

OUTPUT_DIR="${WORKINK_DIR_PREFIX}_anyEpi_outputs"
INPUT_READS_FILE_NAME=$"reads_barcode${BARCODE}"
OUTPUT_FILE="${OUTPUT_DIR}/aligned_${INPUT_READS_FILE_NAME}"

FASTA_FILE_PATH="${OUTPUT_DIR}/fastq/${INPUT_READS_FILE_NAME}.fasta"

#ROOT_DIR_TOOLS="${ROOT_DIR}/tools"
NANOPOLISH_PATH="nanopolish" #"${ROOT_DIR_TOOLS}/nanopolish/nanopolish"
MINIMAP_PATH="minimap2" #${ROOT_DIR_TOOLS}/minimap2/minimap2"
BED_TOOLS_PATH="bedtools"

ROOT_DIR_PYTHON="./"

ANYEPI_ENV_NAME="5gmc_n3_detector"
if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  echo "activate ${ANYEPI_ENV_NAME} conda environment"
  conda activate $ANYEPI_ENV_NAME
else
  echo "${CONDA_DEFAULT_ENV} conda environment is active"
fi

if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  raise error "Error: ${ANYEPI_ENV_NAME} conda environment activation Failed!"
  exit -1
fi

STRAND_SIGN=""
for STRAND in "reverse" "forward"
do

	if [ $STRAND = "forward" ]; then 
		STRAND_SIGN="+"
	else
		STRAND_SIGN="-"	
	fi

	echo  -e "\n\n[-------------]"
	echo "[     <${BARCODE}>     ]"
	echo "[   ${STRAND}   ]"
	echo "[     (${STRAND_SIGN})     ]"
	echo "[-------------]"
	BAM_FILE_PATH="${OUTPUT_FILE}.${STRAND}.sorted.bam"
	OUTPUT_FILE_BASE_PATH="${OUTPUT_FILE}_${STRAND}_${START}_${END}"
	NEW_BAM_FILE_PATH="${OUTPUT_FILE_BASE_PATH}.bam"
	echo -e "\n\nCreate ${STRAND} bam file"
  if [ -f ${NEW_BAM_FILE_PATH} ]; then
		echo -e "[${NEW_BAM_FILE_PATH}] exists"
	else
    samtools view ${OUTPUT_FILE}.${STRAND}.sorted.bam ${CHR}:${START}-${END} -b > ${NEW_BAM_FILE_PATH}
    samtools index  ${NEW_BAM_FILE_PATH}
    samtools view -h ${NEW_BAM_FILE_PATH}  "${OUTPUT_FILE_BASE_PATH}.sam"
  fi

	NEW_FASTA_FILE_PATH="${OUTPUT_FILE_BASE_PATH}.fasta"
	echo -e "\n\nCreate ${STRAND} ${CHR}:${START}-${END} fasta file"
  if [ -f ${NEW_FASTA_FILE_PATH} ]; then
		echo -e "[${NEW_FASTA_FILE_PATH}] exists"
	else
    samtools faidx ${FASTA_FILE_PATH} ${CHR}:${START}-${END} > ${NEW_FASTA_FILE_PATH}
  fi

	EVENTALIGN_FILE_PATH="${OUTPUT_FILE_BASE_PATH}.eventalign.tsv"
	CONSOLIDATED_EVENTS_OUTPUT_FILE="_consolidated_events_${STRAND}_strand.tsv"
	CONSOLIDATED_EVENTALIGN_FILE_PATH="${OUTPUT_FILE_BASE_PATH}${CONSOLIDATED_EVENTS_OUTPUT_FILE}"
	echo -e "\n\nCreate ${STRAND} eventalign file"
	if [ -f ${EVENTALIGN_FILE_PATH} ]; then
		echo -e "[${EVENTALIGN_FILE_PATH}] exists"
	else 
		echo "------"
		echo "${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH}"
		#${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH}
		echo "------"
		echo "${NANOPOLISH_PATH} eventalign --reads ${FASTA_FILE_PATH} --bam ${NEW_BAM_FILE_PATH} --genome ${REF_SEQ} --scale-events --progress --print-read-names -t 16 > ${EVENTALIGN_FILE_PATH}"
		echo "------"
		${NANOPOLISH_PATH} eventalign --reads ${FASTA_FILE_PATH} --bam ${NEW_BAM_FILE_PATH} --genome ${REF_SEQ} --scale-events --progress --print-read-names -t 16 > ${EVENTALIGN_FILE_PATH}
	fi

  echo -e "\n\nCreate ${STRAND} consolidated events eventalign file"
	if [ -f ${CONSOLIDATED_EVENTALIGN_FILE_PATH} ]; then
		echo -e "[${CONSOLIDATED_EVENTALIGN_FILE_PATH}] exists"
	else
		echo "------"
		echo "python ${ROOT_DIR_PYTHON}/nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${EVENTALIGN_FILE_PATH} --ref ${REF_SEQ} -o ${CONSOLIDATED_EVENTALIGN_FILE_PATH} -sp ${START} -ep ${END} -f -s ${STRAND_SIGN} --input-pkl ${PKL_FILE_MODEL}"
		echo "------"
		python ${ROOT_DIR_PYTHON}/nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${EVENTALIGN_FILE_PATH} --ref ${REF_SEQ} -o ${CONSOLIDATED_EVENTALIGN_FILE_PATH} -sp ${START} -ep ${END} -f -s ${STRAND_SIGN} --input-pkl ${PKL_FILE_MODEL}
	fi
done
