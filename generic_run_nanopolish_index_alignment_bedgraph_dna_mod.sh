BARCODE=$1
EXPERIMENT_ROOT_DIR=$2
REF_SEQ=$3
INPUT_READS_FILE_NAME="reads_barcode${BARCODE}"

WORKINK_DIR_PREFIX="${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}"
FAST5_FILES_DIR="${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}/"

OUTPUT_DIR="${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}_epi_basecaller_outputs"

FASTQ_FILES_DIR="${OUTPUT_DIR}/fastq/"
FASTQ_FILE_PATH="${FASTQ_FILES_DIR}/${INPUT_READS_FILE_NAME}.fastq"
FASTA_FILE_PATH="${FASTQ_FILES_DIR}/${INPUT_READS_FILE_NAME}.fasta"

OUTPUT_FILE="${OUTPUT_DIR}/aligned_reads_barcode${BARCODE}"

ROOT_DIR_TOOLS="${ROOT_DIR}/tools"
NANOPOLISH_PATH="nanopolish"
MINIMAP_PATH="minimap2"
BED_TOOLS_PATH=bedtools

ANYEPI_ENV_NAME="smash_seq_env"
if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  echo "conda activate ${ANYEPI_ENV_NAME}"
  conda activate $ANYEPI_ENV_NAME
else
  echo "${CONDA_DEFAULT_ENV} conda environment is active"
fi

echo "Create outputs directory: ${OUTPUT_DIR}"
if [ -d "${OUTPUT_DIR}" ]; then
    echo -e "[${OUTPUT_DIR}] exists"
    #exit -1
else
    mkdir -p ${OUTPUT_DIR}
fi

echo -e "\nPhase 1: Basecall nanopore output fast5 files into fastq files"
if [ -d "${FASTQ_FILES_DIR}" ]; then
    echo -e "[${FASTQ_FILES_DIR}] exists"
else
    if ! [ -d "${FAST5_FILES_DIR}" ]; then
	    echo -e "[${FAST5_FILES_DIR}] does not exists"
      exit -1
	  else
      mkdir -p ${FASTQ_FILES_DIR}
      /opt/ont/ont-guppy/bin/guppy_basecaller --input_path ${FAST5_FILES_DIR} --save_path ${FASTQ_FILES_DIR} \
          --device 'cuda:all:100%' -c dna_r9.4.1_e8.1_fast.cfg # --flowcell FLO-MIN106 --kit SQK_LSK109
    fi
fi

echo -e "\nPhase 2: Create one big fastq file from all the fastq files (the basecalling output)"
if [ -f "${FASTQ_FILE_PATH}" ]; then
    echo -e "[${FASTQ_FILE_PATH}] exists"
else 
    cat ${FASTQ_FILES_DIR}/*/*.fastq >> ${FASTQ_FILE_PATH}
fi

if [ ! -s "${FASTQ_FILE_PATH}" ]; then
        # The file is empty.
	    echo -e "[${FASTQ_FILE_PATH}] is empty"
	    rm -f ${FASTQ_FILE_PATH}
      exit -1
fi

echo -e "\nPhase 2: convert big fastq to fasta format"
if [ -f "${FASTA_FILE_PATH}" ]; then
    echo -e "[${FASTA_FILE_PATH}] exists"
else 
    paste - - - - < ${FASTQ_FILE_PATH} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${FASTA_FILE_PATH}
fi

echo -e "\nPhase 3: Run nanopolish index on the fast5 files directory together with the 'all reads' fasta file"
if [ -f "${FASTA_FILE_PATH}.index" ]; then
    echo -e "[${FASTA_FILE_PATH}.index] exists"
	if [ -f "${FASTA_FILE_PATH}.index.readdb" ]; then
		echo -e "[${FASTA_FILE_PATH}.index.readdb] exists"
	else
		${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH} --verbose --sequencing-summary ${FASTQ_FILES_DIR}/sequencing_summary.txt
	fi
fi

echo -e "\nPhase 4: Align the reads to the reference"
if [ -f "${OUTPUT_FILE}.sam" ]; then
    echo -e "[${OUTPUT_FILE}.sam] exists"
else 
	if [ -f "${OUTPUT_FILE}.bam" ]; then
		echo -e "[${OUTPUT_FILE}.bam] exists"
	else 
		${MINIMAP_PATH} -ax map-ont -t 24 ${REF_SEQ} ${FASTA_FILE_PATH} > ${OUTPUT_FILE}.sam
	fi
fi

echo -e "\nPhase 5: Convert the sam file to a bam file"
if [ -f "${OUTPUT_FILE}.bam" ]; then
    echo -e "[${OUTPUT_FILE}.bam] exists"
else 
	samtools view -S ${OUTPUT_FILE}.sam -b > ${OUTPUT_FILE}.bam
fi

echo -e "\nPhase 5: Sort the bam file"
if [ -f "${OUTPUT_FILE}.sorted.bam" ]; then
    echo -e "[${OUTPUT_FILE}.sorted.bam] exists"
else 
	samtools sort ${OUTPUT_FILE}.bam -o ${OUTPUT_FILE}.sorted.bam
	samtools index ${OUTPUT_FILE}.sorted.bam
fi

echo -e "\nCreate sam file from bam"
if [ -f "${OUTPUT_FILE}.sorted.sam" ]; then
  echo -e "[${OUTPUT_FILE}.sorted.sam] exists"
else
  samtools view -h ${OUTPUT_FILE}.sorted.bam > ${OUTPUT_FILE}.sorted.sam
fi

echo -e "\nCreate strand specific sam files"
if [ -f "${OUTPUT_FILE}.forward.sam" ] && [ -f "${OUTPUT_FILE}.reverse.sam" ]; then
    echo -e "[${OUTPUT_FILE}.forward.sam] exists"
    echo -e "[${OUTPUT_FILE}.reverse.sam] exists"
else
  header=$(cat "${OUTPUT_FILE}.sorted.sam" | grep '^[^@]' --line-number  | head -n 1 | cut -d':' -f1)
  echo "sorted.sam header length=${header}"
  awk -v output_file="${OUTPUT_FILE}" -v header="${header}" '{if(NR<header) {print >> output_file".forward.sam"; print >> output_file".reverse.sam";} if(NR>=header && $2==0) {print $0 >> output_file".forward.sam"} if(NR>=header && $2==16) {print $0 >> output_file".reverse.sam"}}' ${OUTPUT_FILE}.sorted.sam
fi

echo -e "\nCreate forward reads bam file"
if [ -f "${OUTPUT_FILE}.forward.sorted.bam" ]; then
    echo -e "[${OUTPUT_FILE}.forward.sorted.bam] exists"
else
	samtools view -S ${OUTPUT_FILE}.forward.sam -b > ${OUTPUT_FILE}.forward.bam
	samtools sort ${OUTPUT_FILE}.forward.bam -o ${OUTPUT_FILE}.forward.sorted.bam > ${OUTPUT_FILE}.forward.sorted.bam
	samtools index ${OUTPUT_FILE}.forward.sorted.bam
fi

echo -e "\nCreate reverse reads bam file"
if [ -f "${OUTPUT_FILE}.reverse.sorted.bam" ]; then
    echo -e "[${OUTPUT_FILE}.reverse.sorted.bam] exists"
else
	samtools view -S ${OUTPUT_FILE}.reverse.sam -b > ${OUTPUT_FILE}.reverse.bam
	samtools sort ${OUTPUT_FILE}.reverse.bam -o ${OUTPUT_FILE}.reverse.sorted.bam > ${OUTPUT_FILE}.reverse.sorted.bam
	samtools index ${OUTPUT_FILE}.reverse.sorted.bam
fi

echo -e "\nCreate bedgraph files"
if [ -f "${OUTPUT_FILE}.bedgraph" ]; then
    echo -e "[${OUTPUT_FILE}.bedgraph] exists"
else
	${BED_TOOLS_PATH} genomecov -ibam ${OUTPUT_FILE}.sorted.bam -bg > ${OUTPUT_FILE}.bedgraph
	${BED_TOOLS_PATH} genomecov -ibam ${OUTPUT_FILE}.forward.sorted.bam -bg > ${OUTPUT_FILE}.forward.bedgraph
	${BED_TOOLS_PATH} genomecov -ibam ${OUTPUT_FILE}.reverse.sorted.bam -bg > ${OUTPUT_FILE}.reverse.bedgraph
fi