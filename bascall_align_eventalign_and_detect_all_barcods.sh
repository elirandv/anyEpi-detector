#!/bin/bash
EXPERIMENT_ROOT_DIR=$1
BACK_WIND=$2
FRONT_WIN=10
BARCODES=$3
ANYEPI_ENV_NAME="smash_seq_env"
REF_SEQ=/home/nano/Documents/ref_genomes/mouse_c57_bl6/mouse_c57_bl6.sorted.fa #/home/nano/Documents/ref_genomes/e-coli_MG1655/Escherichia_coli_MG1655.fasta
MODEL="./models/random_forest_classifier_V170921.pkl"
modified_sites="/home/nano/Documents/nanopore_piplines/timps_data/test_pipline/barcode02/BL6MBgDNA.CG_5x_only_fully_hmc_in_both_in_genes.bedgraph"
OUTPUT_EVENTALIGN_BASE_PATH=$4

#conda init bash
#exec bash

echo  "conda activate smash_seq_env"
#conda activate $ANYEPI_ENV_NAME

if [ "${CONDA_DEFAULT_ENV}" != "${ANYEPI_ENV_NAME}" ]; then
  raise error "Error: ${ANYEPI_ENV_NAME} conda environment activation Failed!"
  exit -1
else
  echo "${CONDA_DEFAULT_ENV} conda environment is active"
fi

for barcode in $BARCODES
do
    echo "<-----${barcode}----->"   

  echo "generic_run_nanopolish_index_alignment_bedgraph_dna_mod"
	source ./generic_run_nanopolish_index_alignment_bedgraph_dna_mod.sh ${barcode} ${EXPERIMENT_ROOT_DIR} ${REF_SEQ}

  while read -r CHR b c OTHER
  do
    reverse="reverse"
    forward="forward"
    STRAND_2=$(echo $OTHER | sed 's/ .*//')
    if [[ "$STRAND_2" != "$forward" && "$STRAND_2" != "$reverse" ]];
    then
      echo "No strand was mentioned, Running ${CHR} ${b}-${c} as forward & reverse location"
    else
      reverse=""
      forward=$STRAND_2
    fi
    for STRAND in $forward $reverse
    do
      echo $CHR
#    done
#  done < "${EXPERIMENT_ROOT_DIR}/barcode03/modified_sites_F_R-5+5_kmer_lex_like_kmer_2881672_2881700.bedgraph" # "${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}/modified_sites.bedgraph"
#done
#exit
      back_window=$4
      START=$((b-BACK_WIND));
      END=$((b+FRONT_WIN));
      echo  -e "\n\n[--------------------]"
      echo "[        <${BARCODE}>       "
      echo "[       ${CHR}       "
      echo "[      ${STRAND}      "
      echo "[  ${START}-${END}"
      echo "[--------------------]"

      echo "create_bam_and_eventalign_files_for_specific_region_dna_mod"
      source ./create_bam_and_eventalign_files_for_specific_region_dna_mod.sh ${EXPERIMENT_ROOT_DIR} ${REF_SEQ} ${MODEL} ${barcode} ${CHR} ${START} ${END} ${STRAND}  $OUTPUT_EVENTALIGN_BASE_PATH
    #  echo "cut_eventalign_and_run_detector_for_bed_selected_ereas"
    #  source ./cut_eventalign_and_run_detector_for_bed_selected_ereas.sh ${modified_sites} ${barcode}
    done
  done < "${EXPERIMENT_ROOT_DIR}/barcode03/U_BED_CG_RF" # "${EXPERIMENT_ROOT_DIR}/barcode${BARCODE}/modified_sites.bedgraph"
done
