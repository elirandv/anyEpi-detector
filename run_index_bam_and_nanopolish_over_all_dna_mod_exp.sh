ROOT_DIR="/ebensteinlab"
SCRIPT_DIR="${ROOT_DIR}/projects/DNA_modifications_030620/scripts/bash"
PROJECT_ROOT_DIR="${ROOT_DIR}/projects/DNA_modifications_030620/DNA_modifications_030620/20200603_1039_MN18252_FAN25046_572476a1"

ROOT_DIR_TOOLS="${ROOT_DIR}/tools"
NANOPOLISH_PATH="${ROOT_DIR_TOOLS}/nanopolish/nanopolish"

module load python/python-3.6.7

for barcode in 1 2 3 4 5
do 	
    echo "<-----${barcode}----->"   
	
	ROOT_DIR_RUN_FILES="${PROJECT_ROOT_DIR}/barcode0${barcode}"
	FAST5_FILES_DIR="${ROOT_DIR_RUN_FILES}_fast5"
	FASTQ_FILES_DIR="${ROOT_DIR_RUN_FILES}_guppy4_fastq"

	OUTPUT_DIR="${ROOT_DIR_RUN_FILES}_outputs"
	INPUT_READS_FILE_NAME="reads_barcode0${barcode}"
	FASTA_FILE_PATH="${OUTPUT_DIR}/${INPUT_READS_FILE_NAME}.fasta"
	
	echo "Index the fasta reads file"
	if [ -f "${FASTA_FILE_PATH}.index" ]; then
		echo -e "[${FASTA_FILE_PATH}.index] exists"
		if [ -f "${FASTA_FILE_PATH}.index.readdb" ]; then
			echo -e "[${FASTA_FILE_PATH}.index.readdb] exists"
		else 
			${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH} --verbose 
		fi
	else
		${NANOPOLISH_PATH} index -d ${FAST5_FILES_DIR} ${FASTA_FILE_PATH} --verbose 
	fi
	
	${SCRIPT_DIR}/create_bam_and_eventalign_files_for_specific_region_dna_mod.sh ${barcode} 10000 11060
	
	${SCRIPT_DIR}/create_bam_and_eventalign_files_for_specific_region_dna_mod.sh ${barcode} 39610 42430
done
