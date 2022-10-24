#!/usr/bin/mawk -f 
BEGIN{OFS="\t";
      removeSelfHit = 1;
      tp_fam = 0; tp_sfam = 0; tp_fold = 0; fp = 0;	
      print "PREC_FAM","PREC_SFAM","PREC_FOLD","RECALL_FAM","RECALL_SFAM","RECALL_FOLD";
}
FNR==NR{
	id2fam[$1]=$2; 
        famCnt[$2]++; 
        gsub(/\.[0-9]+$/,"",$2); 
        id2sfam[$1]=$2; 
        sfamCnt[$2]++; 
        gsub(/\.[0-9]+$/,"",$2); 
        id2fold[$1]=$2; 
        foldCnt[$2]++; 
        next }
removeSelfHit == 1 { removeSelfHit = 0;
                queries = 0;
	        for(fam in famCnt){
		    famCnt[fam] = famCnt[fam] - 1;
		}
	        for(sfam in sfamCnt){
		    sfamCnt[sfam] = sfamCnt[sfam] - 1;
		}
		for(fold in foldCnt){
		    foldCnt[fold] = foldCnt[fold] - 1;
		}
		for(id in id2fam){
		  count =  (famCnt[id2fam[id]] != 0 && sfamCnt[id2sfam[id]] - famCnt[id2fam[id]] != 0 && foldCnt[id2fold[id]] - sfamCnt[id2sfam[id]] != 0);
		  queries = queries + count;
	        } 
}	
!($1 in id2fam) {next}	
!($2 in id2fam) {next}
$1 == $2 {next} # skip self hit

NR %1000 == 0{ print tp_fam/(tp_fam+fp), tp_sfam/(tp_sfam+fp), tp_fold/(tp_fold+fp), tp_fam / queries, tp_sfam / queries, tp_fold / queries, tp_fam; }										       
id2fold[$1] != id2fold[$2] {norm=(foldCnt[id2fold[$1]] - sfamCnt[id2sfam[$1]]); norm = (norm > famCnt[id2fam[$1]]) ? norm : famCnt[id2fam[$1]]; norm = (norm == 0) ? 1 : norm; fp = fp + (1 / norm); next }
(famCnt[id2fam[$1]] == 0 || sfamCnt[id2sfam[$1]] - famCnt[id2fam[$1]] == 0 || foldCnt[id2fold[$1]] - sfamCnt[id2sfam[$1]] == 0){ next } 

id2fam[$1] == id2fam[$2] { norm=famCnt[id2fam[$1]];
	                   tp_fam = tp_fam + (1 / norm);
			   next }
id2fam[$1] != id2fam[$2] && id2sfam[$1] == id2sfam[$2]  { norm=(sfamCnt[id2sfam[$1]] - famCnt[id2fam[$1]]);
	                                                  tp_sfam = tp_sfam + (1 / norm);
							  next }
id2fam[$1] != id2fam[$2] && id2sfam[$1] != id2sfam[$2] && id2fold[$1] == id2fold[$2] { norm=(foldCnt[id2fold[$1]] - sfamCnt[id2sfam[$1]]); 
	                                                                               tp_fold = tp_fold + (1 / norm); 
										       next }

END{
       print tp_fam/(tp_fam+fp), tp_sfam/(tp_sfam+fp), tp_fold/(tp_fold+fp), tp_fam / queries, tp_sfam / queries, tp_fold / queries, tp_fam; 									       
       #print tp_fam, tp_sfam, tp_fold, fp;	
}
