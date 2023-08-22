back_window=$2
front_window=10
while read -r a b c d e
do
  sp=$((b-back_window));
  ep=$((b+front_window));
  echo "<----- -sp ${sp} -ep ${ep} ----->";

  STRAND="forward"
  BARCODES="09,10,11,12"
  WORK_DIR="./220299_5hmc_training/"
  COUNT_EVENTS_AND_NNN_PLOT="${WORK_DIR}_${sp}_${ep}_events_count.png"

  echo -e "\n\n plot events and nnn count and ratio from ${sp} to ${ep}"
  if [ -f ${COUNT_EVENTS_AND_NNN_PLOT} ]; then
    echo -e "[${COUNT_EVENTS_AND_NNN_PLOT}] exists"
  else
    echo "------"
    echo "python ./k-mer_plot_events_and_nnn_count.py -o ${WORK_DIR} -b ${BARCODES} -s ${STRAND} -seq ${d} -i ${sp}_${ep} -w ${back_window}"
    echo "------"
    python ./k-mer_plot_events_and_nnn_count.py -o ${WORK_DIR} -b ${BARCODES} -s ${STRAND} -seq ${d} -i ${sp}_${ep} -w ${back_window}
  fi
done < $1