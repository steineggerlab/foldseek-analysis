#!/usr/bin/mawk -f 
BEGIN{OFS="\t";
      print "NAME","SCOP","FAM","SFAM","FOLD","FP","FAMCNT","SFAMCNT","FOLDCNT";
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
foundFp[$1] < 1 && id2fold[$1] != id2fold[$2] {foundFp[$1]++; next}
foundFp[$1] < 1 && id2fam[$1] == id2fam[$2] { foundFam[$1]++; next }
foundFp[$1] < 1 && id2fam[$1] != id2fam[$2] && id2sfam[$1] == id2sfam[$2] { foundSFam[$1]++; next }
foundFp[$1] < 1 && id2fam[$1] != id2fam[$2] && id2sfam[$1] != id2sfam[$2] && id2fold[$1] == id2fold[$2] { foundFold[$1]++; next }
END{ 
   for(i in id2fam){
      if(id2fam[i] != "" && sfamCnt[id2sfam[i]] - famCnt[id2fam[i]] > 0 && foldCnt[id2fold[i]] - sfamCnt[id2sfam[i]] > 0){	   
        famVal=foundFam[i]/famCnt[id2fam[i]];
        sfamVal=foundSFam[i]/(sfamCnt[id2sfam[i]] - famCnt[id2fam[i]]);
        foldVal=foundFold[i]/(foldCnt[id2fold[i]] - sfamCnt[id2sfam[i]]);
	fpCnt = (foundFp[i] == "") ? 0 : foundFp[i]; 
        print i,id2fam[i],famVal,sfamVal,foldVal,fpCnt,famCnt[id2fam[i]],sfamCnt[id2sfam[i]],foldCnt[id2fold[i]];
      }
   }  
}
