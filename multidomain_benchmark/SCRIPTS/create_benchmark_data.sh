# SPICi clustering
## convert alphafold database to fasta (foldseek Version: fdc503ea7fb6c34991146ccf9bdec07125f1b6ba)
foldseek convert2fasta alphafolddb_v1 db.fasta
## use BLAST (makeblastdb: 2.5.0+)
makeblastdb -in db.fasta -out blast_db -dbtype prot
# Protein-Protein BLAST 2.5.0+
blastp -db blast_db -query db.fasta -out blast_results.out
awk '$11 < 0.001 {print $1 "\t" $2 "\t" $12}' blast_results.out > blast_lowerE_results.out
## use spici (use normalised score) (Extremely fast graph clustering algorithm. (Peng Jiang 2009))
# downloaded on the  Apr 11 2022 from https://compbio.cs.princeton.edu/spici/
spici -i blast_lowerE_results.out -o spici_lowerE.txt
# pick the longest sequence for each cluster
python3 pick_longest_spici.py spici_lowerE.txt > spici_longest.txt
# use spici results to make new target
for i in $(cat spici_longest.txt); do cp ./pdbfiles/$i spici_target/$i; done
# use spici results to make new query
## pick 100 at random to be the queries
cat spici_longest.txt | shuf | head -n 100 > spici_query_100_random.txt
## cope them to query folder
for i in $(cat spici_query_100_random.txt); do cp spici_target/$i spici_query/; done
# remove queries from target
for i in *; do (rm ../spici_target/$i); done
