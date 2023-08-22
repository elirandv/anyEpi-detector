EXPERIMENT_ROOT_DIR=$1
REF_SEQ=$2
PKL_FILE_MODEL=$3
BARCODE=$4
CHR=$5
START=$6
END=$7
INPUT_READS_FILE_NAME="reads_barcode${BARCODE}"
STRAND=$8
OUTPUT_EVENTALIGN_BASE_PATH=$9

WORKINK_DIR_PREFIX="${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}"
FAST5_FILES_DIR="${WORKINK_DIR_PREFIX}/"

OUTPUT_DIR="${WORKINK_DIR_PREFIX}_epi_basecaller_outputs"
OUTPUT_FILE="${OUTPUT_DIR}/aligned_${INPUT_READS_FILE_NAME}"

FASTA_FILE_PATH="${OUTPUT_DIR}/fastq/${INPUT_READS_FILE_NAME}.fasta"

#ROOT_DIR_TOOLS="${ROOT_DIR}/tools"
NANOPOLISH_PATH="/home/nano/miniconda3/envs/smash_seq_env/bin/nanopolish/nanopolish" # "nanopolish"
MINIMAP_PATH="minimap2" #${ROOT_DIR_TOOLS}/minimap2/minimap2"
BED_TOOLS_PATH="bedtools"

ROOT_DIR_PYTHON="./"

ANYEPI_ENV_NAME="smash_seq_env"
if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  echo "conda activate ${ANYEPI_ENV_NAME}"
  conda activate $ANYEPI_ENV_NAME
else
  echo "${CONDA_DEFAULT_ENV} conda environment is active"
fi

if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  echo "Error: ${ANYEPI_ENV_NAME} conda environment activation Failed!"
fi

STRAND_SIGN=""
if [ $STRAND = "forward" ]; then
  STRAND_SIGN="+"
else
  STRAND_SIGN="-"
fi

BAM_FILE_PATH="${OUTPUT_FILE}.${STRAND}.sorted.bam"

OUTPUT_FILE_BASE_PATH="${OUTPUT_EVENTALIGN_BASE_PATH}/barcode${BARCODE}_epi_basecaller_outputs/eventalign_${STRAND}/"
echo ${OUTPUT_EVENTALIGN_BASE_PATH}
#	exit
FILE_NAME_BASE="${OUTPUT_FILE_BASE_PATH}/${START}_${END}"
#FILE_NAME_BASE="${OUTPUT_EVENTALIGN_BASE_PATH}aligned_${INPUT_READS_FILE_NAME}_${START}_${END}_${STRAND}"
NEW_BAM_FILE_PATH="${FILE_NAME_BASE}.bam"
echo -e "\n\nCreate ${STRAND} bam file"
if [ -f ${NEW_BAM_FILE_PATH} ]; then
  echo -e "[${NEW_BAM_FILE_PATH}] exists"
else
  samtools view ${OUTPUT_FILE}.${STRAND}.sorted.bam ${CHR}:${START}-${END} -b > ${NEW_BAM_FILE_PATH}
  samtools index  ${NEW_BAM_FILE_PATH}
  samtools view -h ${NEW_BAM_FILE_PATH}  "${FILE_NAME_BASE}.sam"
fi

# try to take smaller fasta
#	NEW_FASTA_FILE_PATH="${FILE_NAME_BASE}.fasta"
#	echo -e "\n\nCreate ${STRAND} ${CHR}:${START}-${END} fasta file"
#  if [ -f ${NEW_FASTA_FILE_PATH} ]; then
#		echo -e "[${NEW_FASTA_FILE_PATH}] exists"
#	else
#    samtools faidx ${FASTA_FILE_PATH} ${CHR}:${START}-${END} -o ${NEW_FASTA_FILE_PATH}
#  fi

if [ -f "${FASTA_FILE_PATH}.index.readdb" ]; then
  echo -e "[${FASTA_FILE_PATH}.index.readdb] exists"
else
  echo "------"
  echo "${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH}"
  ${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH}
  echo "------"
fi

echo "Create outputs directory: ${OUTPUT_FILE_BASE_PATH}"
if [ -d "${OUTPUT_FILE_BASE_PATH}" ]; then
    echo -e "[${OUTPUT_FILE_BASE_PATH}] exists"
else
    mkdir -p ${OUTPUT_FILE_BASE_PATH}
fi

EVENTALIGN_FILE_PATH="${FILE_NAME_BASE}.eventalign.tsv"
CONSOLIDATED_EVENTALIGN_FILE_PATH="${FILE_NAME_BASE}_consolidated_events.tsv"

echo -e "\n\nCreate ${STRAND} eventalign file"
if [ -f "${FILE_NAME_BASE}${START}_${END}.eventalign.tsv" ]; then
  echo -e "["${FILE_NAME_BASE}${START}_${END}.eventalign.tsv"] exists"
else
  echo "------"
  HDF5_PLUGIN_PATH=$(whereis hdf5 | awk '{print($2)}')
  export HDF5_PLUGIN_PATH=$HDF5_PLUGIN_PATH/lib/plugin/
  TSV_PATH="${OUTPUT_EVENTALIGN_FILES}_${STRAND}"
  echo "${NANOPOLISH_PATH} eventalign --reads ${FASTA_FILE_PATH} --bam ${NEW_BAM_FILE_PATH} --genome ${REF_SEQ} --scale-events --progress --print-read-names -t 16 --samples > ${OUTPUT_FILE_BASE_PATH}temp.eventalign.tsv"
  ${NANOPOLISH_PATH} eventalign --reads ${FASTA_FILE_PATH} --bam ${NEW_BAM_FILE_PATH} --genome ${REF_SEQ} --scale-events --progress --print-read-names -t 16 --samples > "${OUTPUT_FILE_BASE_PATH}temp.eventalign.tsv"
  awk '{if ($2>'${START}' && $2<'${END}') print}' "${OUTPUT_FILE_BASE_PATH}temp.eventalign.tsv" > "${FILE_NAME_BASE}.eventalign.tsv"
  echo "------"
fi
rm $NEW_BAM_FILE_PATH*
rm  ${OUTPUT_FILE_BASE_PATH}temp.eventalign.tsv
#  echo -e "\n\nCreate ${STRAND} consolidated events eventalign file"
#	if [ -f ${CONSOLIDATED_EVENTALIGN_FILE_PATH} ]; then
#		echo -e "[${CONSOLIDATED_EVENTALIGN_FILE_PATH}] exists"
#	else
#		echo "------"
#		echo "python ${ROOT_DIR_PYTHON}/nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${EVENTALIGN_FILE_PATH} --ref ${REF_SEQ} -o ${CONSOLIDATED_EVENTALIGN_FILE_PATH} -sp ${START} -ep ${END} -f -s ${STRAND_SIGN} --input-pkl ${PKL_FILE_MODEL}"
#		echo "------"
#		python ${ROOT_DIR_PYTHON}/nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${EVENTALIGN_FILE_PATH} --ref ${REF_SEQ} -o ${CONSOLIDATED_EVENTALIGN_FILE_PATH} -sp ${START} -ep ${END} -f -s ${STRAND_SIGN} --input-pkl ${PKL_FILE_MODEL}
#	fi
echo "done"
