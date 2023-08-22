BED="/home/nano/Documents/nanopore_piplines/timps_data/test_pipline/barcode03/modified_sites_0.7.bedgraph"
sort $BED | uniq > U_BED

bedtools getfasta -fi '/home/nano/Documents/ref_genomes/mouse_c57_bl6/mouse_c57_bl6.sorted.fa' -bed $U_BED -tab -bedOut -s  -name | awk '{print toupper($5)}' | sort | uniq -c


bedtools getfasta -fi '/home/nano/Documents/ref_genomes/mouse_c57_bl6/mouse_c57_bl6.sorted.fa' -bed $U_BED -tab -bedOut -s  -name | awk 'BEGIN{ OFS="\t"} {if (toupper($5)=="G") print $1,$2,$3,"reverse",$4,"-",toupper($5); else print $1,$2,$3,"forward",$4,"+",toupper($5)}' > U_BED_CG

awk 'BEGIN{ OFS="\t"} {print $1,$2-5,$3+5,$4,$5,$6;}' $U_BED_CG  > U_BED_CG_RF

bedtools getfasta -fi '/home/nano/Documents/ref_genomes/mouse_c57_bl6/mouse_c57_bl6.sorted.fa'   -tab -bedOut -s  -name  -bed $U_BED_CG_RF > $U_BED_CG_RF_slope_5

cat $U_BED_CG_RF_slope_5 | awk 'BEGIN{ OFS="\t"} {split(toupper($7),a,""); print $1, $2, $3, $4, $5, $6, $7, a[5]a[8]a[4]a[9]a[3]a[10]a[2]a[11]a[1]}' | sort -k8  > U_BED_CG_RF_slope_5_kmer_lex

cat $U_BED_CG_RF_slope_5_kmer_lex '/home/nano/Documents/nanopore_piplines/timps_data/test_pipline/barcode03/sinthetic_25_kmers_lex.bedgraph' | sort -k8 | column -t > /home/nano/Documents/nanopore_piplines/timps_data/test_pipline/barcode03/modified_sites_0.7_RF_slope_5_kmer_lex_with_sinthetic_25.bedgraph







