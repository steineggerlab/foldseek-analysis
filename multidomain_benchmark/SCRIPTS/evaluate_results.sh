#!/bin/sh -e
# evaluate multidomain benchmark results
TP=$1
FP=$2
LDDT=$3
HITS=$4
NAME=$5

# Coverage
## Foldseek
python3 lddt_20_plddt_avelddt_lddtcutoff.py --alignment foldseek_multidomain_results --query spici_query/ --target spici_target/--tp-threshold $TP --fp-threshold $FP --lddt-cutoff $LDDT > "foldseek_lddt_cutoff_"$NAME".txt"
# MMseqs2
python3 lddt_20_plddt_avelddt_lddtcutoff.py --alignment  mmseqs_spici_results.m8 --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP --lddt-cutoff $LDDT > "mmseqs_lddt_cutoff_"$NAME".txt"
## TM-align
python3 lddt_20_plddt_avelddt_lddtcutoff.py --alignment tmaln_spici.100.out --query spici_query/ --target spici_target/--tp-threshold $TP --fp-threshold $FP --lddt-cutoff $LDDT > "tmalign_lddt_cutoff_"$NAME".txt"
## DALi
python3 lddt_20_plddt_avelddt_lddtcutoff.py --alignment dali_spici.100.pdb.out --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP --lddt-cutoff $LDDT > "dali_lddt_cutoff_"$NAME".txt"
## CLE-SW
python3 lddt_20_plddt_avelddt_lddtcutoff.py --alignment clesw.spici.100.out --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP --lddt-cutoff $LDDT > "clesw_lddt_cutoff_"$NAME".txt"


# Alignment quality
"TOOL, SENSITIVITY, SELECTIVITY"  >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
## Foldseek
"Foldseek," >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
python3 lddt_20_quali_plddt70.py --alignment foldseek_multidomain_results  --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP > "ALI_QUALI_PLDDT70_FS_"$NAME".txt"
python3 parse_ali_quality.py  "ALI_QUALI_PLDDT70_FS_"$NAME".txt" $TP $HITS >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
## MMseqs2
"MMseqs2," >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
python3 lddt_20_quali_plddt70.py --alignment  mmseqs_spici_results.m8 --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP >  "ALI_QUALI_PLDDT70_MM_"$NAME".txt"
python3 parse_ali_quality.py "ALI_QUALI_PLDDT70_MM_"$NAME".txt" $TP $HITS >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
## TM-align
"TMalign," >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
python3 lddt_20_quali_plddt70.py --alignment tmaln_spici.100.out  --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP > "ALI_QUALI_PLDDT70_TMA_"$NAME".txt"
python3 parse_ali_quality.py "ALI_QUALI_PLDDT70_TMA_"$NAME".txt" $TP $HITS >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
## DALi
"DALI," >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
python3 lddt_20_quali_plddt70.py --alignment dali_spici.100.pdb.out  --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP > "ALI_QUALI_PLDDT70_DALI_"$NAME".txt"
python3 parse_ali_quality.py "ALI_QUALI_PLDDT70_DALI_"$NAME".txt" $TP $HITS >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
## CLE-SW
"CLE-SW," >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
python3 lddt_20_quali_plddt70.py --alignment clesw.spici.100.out  --query spici_query/ --target spici_target/ --tp-threshold $TP --fp-threshold $FP > "ALI_QUALI_PLDDT70_CLE_"$NAME".txt"
python3 parse_ali_quality.py "ALI_QUALI_PLDDT70_CLE_"$NAME".txt" $TP $HITS >> "sens_selec_"$HITS"hit_"$TP"TP.csv"
