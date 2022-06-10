#AlphaFold distribution of first 1/5/20 hits (set hit paraemter accordingly)
## Foldseek
python3 lddt_dis_AF.py --alignment foldseek_multidomain_results  --query spici_query/ --target spici_target/ --hits 5 > fs_af_quali.txt
python3 parse-lddt-distr.py fs_af_quali.txt > AF_FS.txt
## MMseqs2
python3 lddt_dis_AF.py --alignment  mmseqs_spici_results.m8 --query spici_query/ --target spici_target/ --hits 5 > mm_af_quali.txt
python3 parse-lddt-distr.py mm_af_quali.txt > AF_MM.txt
## TM-align
python3 lddt_dis_AF.py --alignment tmaln_spici.100.out  --query spici_query/ --target spici_target/ --hits 5 > tma_af_quali.txt
python3 parse-lddt-distr.py tma_af_quali.txt > AF_TMA.txt
## DALi
python3 lddt_dis_AF.py --alignment dali_spici.100.pdb.out  --query spici_query/ --target spici_target/ --hits 5 > dali_af_quali.txt
python3 parse-lddt-distr.py dali_af_quali.txt > AF_DALI.txt
## CLE-SW
python3 lddt_dis_AF.py --alignment clesw.spici.100.out  --query spici_query/ --target spici_target/ --hits 5 > cle_af_quali.txt
python3 parse-lddt-distr.py cle_af_quali.txt > AF_CLE.txt

# SCOPe distribution
## Foldseek
# TODO STEPH: can you fill in and upload the relevant files?
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > fs_fp_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > fs_tp_fam_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > fs_tp_sfam_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > fs_tp_fold_scop_dist.txt
python3 parse-lddt-distr.py fs_fp_scop_dist.txt > FP_FS_AVERAGE.txt
python3 parse-lddt-distr.py fs_tp_fam_scop_dist.txt > TP_FS_FAM_AVERAGE.txt
python3 parse-lddt-distr.py fs_tp_sfam_scop_dist.txt > TP_FS_SFAM_AVERAGE.txt
python3 parse-lddt-distr.py fs_tp_fold_scop_dist.txt > TP_FS_FOLD_AVERAGE.txt

## TM-align
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > tma_fp_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > tma_tp_fam_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > tma_tp_sfam_scop_dist.txt
python3 lddt_dis.py --alignment --query scop-pdb/ --target scop-pdb/ > tma_tp_fold_scop_dist.txt
python3 parse-lddt-distr.py tma_fp_scop_dist.txt > FP_TMA_AVERAGE.txt
python3 parse-lddt-distr.py tma_tp_fam_scop_dist.txt > TP_TMA_FAM_AVERAGE.txt
python3 parse-lddt-distr.py tma_tp_sfam_scop_dist.txt > TP_TMA_SFAM_AVERAGE.txt
python3 parse-lddt-distr.py tma_tp_fold_scop_dist.txt > TP_TMA_FOLD_AVERAGE.txt

## DALi
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > dali_fp_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > dali_tp_fam_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > dali_tp_sfam_scop_dist.txt
python3 lddt_dis.py --alignment  --query scop-pdb/ --target scop-pdb/ > dali_tp_fold_scop_dist.txt
python3 parse-lddt-distr.py dali_fp_scop_dist.txt > FP_DALI_AVERAGE.txt
python3 parse-lddt-distr.py dali_tp_fam_scop_dist.txt > TP_DALI_FAM_AVERAGE.txt
python3 parse-lddt-distr.py dali_tp_sfam_scop_dist.txt > TP_DALI_SFAM_AVERAGE.txt
python3 parse-lddt-distr.py dali_tp_fold_scop_dist.txt > TP_DALI_FOLD_AVERAGE.txt

## CLE-SW
python3 lddt_dis.py --alignment rand_clesw.fp.10000.list --query scop-pdb/ --target scop-pdb/ > cle_fp_scop_dist.txt
python3 lddt_dis.py --alignment rand_clesw_tp_fam.10000.list --query scop-pdb/ --target scop-pdb/ > cle_tp_fam_scop_dist.txt
python3 lddt_dis.py --alignment rand_clesw_tp_sfam.10000.list --query scop-pdb/ --target scop-pdb/ > cle_tp_sfam_scop_dist.txt
python3 lddt_dis.py --alignment rand_clesw_tp_fold.10000.list --query scop-pdb/ --target scop-pdb/ > cle_tp_fold_scop_dist.txt
python3 parse-lddt-distr.py cle_fp_scop_dist.txt > FP_CLESW_AVERAGE.txt
python3 parse-lddt-distr.py cle_tp_fam_scop_dist.txt > TP_CLESW_FAM_AVERAGE.txt
python3 parse-lddt-distr.py cle_tp_sfam_scop_dist.txt > TP_CLESW_SFAM_AVERAGE.txt
python3 parse-lddt-distr.py cle_tp_fold_scop_dist.txt > TP_CLESW_FOLD_AVERAGE.txt
