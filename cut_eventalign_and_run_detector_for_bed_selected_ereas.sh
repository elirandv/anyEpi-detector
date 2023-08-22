OUTPUT_EVENTALIGN_FILES=$1
BARCODE=$3

while read -r a b c d e
do
  back_window=$4
  sp=$((b-back_window));
  ep=$((b+10));
  echo "<----- -sp ${sp} -ep ${ep} ----->";

  STRAND_SIGN=""
  for STRAND in "forward" #"reverse"
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

    TSV_PATH="${OUTPUT_EVENTALIGN_FILES}_${STRAND}"

    echo "Create outputs directory: ${TSV_PATH}/"
    if [ -d "${TSV_PATH}/" ]; then
        echo -e "[${TSV_PATH}/] exists"
    else
        mkdir -p ${TSV_PATH}/
    fi

#from full tsv
    tsv="${TSV_PATH}.eventalign.tsv"
#from already cut tsv
#    f=$((b-25));
#    tsv="${TSV_PATH}_${f}_${ep}.eventalign.tsv"

    new_tsv="${TSV_PATH}/${sp}_${ep}.eventalign.tsv"
    echo -e "\n\nCut eventalign.tsv from ${sp} to ${ep}"
    if [ -f ${new_tsv} ]; then
      echo -e "[${new_tsv}] exists"
    else
      echo "------"
      echo 'awk {if ($2>'${sp}' && $2<'${ep}') print}' "${tsv} > ${new_tsv}"
      echo "------"
      awk '{if ($2>'${sp}' && $2<'${ep}') print}' ${tsv} > ${new_tsv}
    fi

#    COUNT_EVENTS_AND_NNN_TSV="${TSV_PATH}/_${sp}_${ep}_events_count.tsv"
#    echo -e "\n\n Count events and nnn ratio from ${sp} to ${ep}"
#    if [ -f ${COUNT_EVENTS_AND_NNN_TSV} ]; then
#      echo -e "[${COUNT_EVENTS_AND_NNN_TSV}] exists"
#    else
#      echo "------"
#      echo "awk  {if ($10=="NNNNNN") print $2, 0; else print $2, 1} ${new_tsv} | sort | uniq -c"
#      echo "------"
#      awk ' {if ($10=="NNNNNN") print $2, 0; else print $2, 1} ' ${new_tsv} | sort | uniq -c | awk 'NR%2{printf "%s ",$0;next;}1' | awk -v VAR=$back_window 'BEGIN {loc=-VAR+1; print "mod_dist","event_k-mer_first_loc", "event_count", "NNN_count", "NNN_ratio"} {print loc++, $2, $1+$4, $1, $1/($1+$4)}'  > ${COUNT_EVENTS_AND_NNN_TSV}
#    fi

    COUNT_EVENTS_AND_NNN_TSV="${TSV_PATH}/${sp}_${ep}_events_count_norm.tsv"
    echo -e "\n\n Count events and nnn ratio from ${sp} to ${ep}"
    if [ -f ${COUNT_EVENTS_AND_NNN_TSV} ]; then
      echo -e "[${COUNT_EVENTS_AND_NNN_TSV}] exists"
    else
      echo "------"
      echo "awk  {if ($10=="NNNNNN") print $2, 0; else print $2, 1} ${new_tsv} | sort | uniq -c"
      echo "------"
      awk 'BEGIN {for (i='${sp}';i<'${ep}';i++) print i,0,0"\n"i,1,0 } {if ($10=="NNNNNN") print $2, 0, $9; else print $2, 1, $9} ' ${new_tsv} | awk '{ a[$1,$2] += $3 }
  END {
    for (i in a) {
      printf "%-15s\t%s\n", a[i], i | "sort -k 2";
    }
  }' | awk 'NR%2{printf "%s ",$0;next;}1' | awk '{sum+=$1+$3} {print} END {print sum}' | awk -v VAR=$back_window 'BEGIN {loc=-VAR+1; print "mod_dist","event_k-mer_first_loc", "event_count", "event_norm", "NNN_count", "NNN_ratio"} {print loc++, $2, $1+$3, $1, $1/($1+$3)}' >  "${TSV_PATH}_temp.tsv"
      sum_events=$(tail -1 "${TSV_PATH}_temp.tsv" | awk '{print $2}')
      awk -v VAR=$sum_events 'NR==1 {print} {if ($1<10) print $1, substr($2, 1, length($2)-2), $3, $3/VAR, $4, $5}' "${TSV_PATH}_temp.tsv"  > ${COUNT_EVENTS_AND_NNN_TSV}
    fi

    #COUNT_EVENTS_AND_NNN_PLOT="${TSV_PATH}_${sp}_${ep}_events_count.png"
    #echo -e "\n\n plot events and nnn count and ratio from ${sp} to ${ep}"
    #if [ -f ${COUNT_EVENTS_AND_NNN_PLOT} ]; then
    #  echo -e "[${COUNT_EVENTS_AND_NNN_PLOT}] exists"
    #else
    #  echo "------"
    #  echo "python ./plot_events_and_nnn_count.py -i ${COUNT_EVENTS_AND_NNN_TSV} -o ${COUNT_EVENTS_AND_NNN_PLOT} -b ${BARCODE} -s ${STRAND}  -l ${sp}_${ep} -seq $(d) -w ${back_window}"
    #  echo "------"
    #  python ./plot_events_and_nnn_count.py -i ${COUNT_EVENTS_AND_NNN_TSV} -o ${COUNT_EVENTS_AND_NNN_PLOT} -b ${BARCODE} -s ${STRAND} -l ${sp}_${ep} -seq $(d) -w ${back_window}
    #fi

#    out_tsv="${TSV_PATH}_${sp}_${ep}_consolidated_events.tsv"
#    echo -e "\n\nCreate ${STRAND} consolidated events eventalign file"
#    if [ -f ${out_tsv} ]; then
#      echo -e "[${out_tsv}] exists"
#    else
#      echo "------"
#      echo "python ./nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${new_tsv} --ref /home/nano/Documents/ref_genomes/e-coli_MG1655/Escherichia_coli_MG1655.fasta -o ${out_tsv} -sp ${sp} -ep ${ep} -f -s ${STRAND_SIGN} --input-pkl ./models/random_forest_classifier_V170921.pkl"
#      echo "------"
#      python ./nanopore_sequencing_5gmc_n3_detector_v1.minion_dna_r941.py extract-features-and-detect -i ${new_tsv} --ref /home/nano/Documents/ref_genomes/e-coli_MG1655/Escherichia_coli_MG1655.fasta -o ${out_tsv} -sp ${sp} -ep ${ep} -f -s ${STRAND_SIGN} --input-pkl ./models/random_forest_classifier_V170921.pkl
#    fi
  done
done < $2 # ./220299_5hmc_training/modified_sites_both_strands_and_11mer.bedgraph
