#!/bin/bash -xe

# Build Foldseek
git clone git@github.com:steineggerlab/foldseek.git && cd foldseek
git checkout fdc503ea7fb6c34991146ccf9bdec07125f1b6ba
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
patch -u ../src/strucclustutils/tmalign.cpp -i ../../foldseek_measurment.patch
make -j
make install
cd ../..

# Fetch some pdbs
mkdir -p pdbs

PDBS=$(find ../../../data/structur_align/pdbs/ -name 'd*' | shuf | head -n 1000)
for i in $PDBS; do
  cp $i pdbs/$(basename $i)
done

# Run Foldseek
mkdir -p tmp

foldseek/build/src/foldseek createdb pdbs tmp/queryDB --threads 1

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "${TDB}.index" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
    awk -v size=$INDEX_SIZE '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

rm tmp/fakeprefDB*  # TODO
fake_pref tmp/queryDB tmp/queryDB tmp/fakeprefDB
ln -sf queryDB.index tmp/fakeprefDB

foldseek/build/src/foldseek tmalign tmp/queryDB tmp/queryDB tmp/fakeprefDB tmp/tmalignRes \
  --threads 1 --tmscore-threshold 0 --tmalign-fast 1 \
  > tmp/foldseek_output.txt

awk 'FNR==NR{id2sid[$1]=$2;next}
     /myduration/{print id2sid[$6],id2sid[$7],$2,$3,$4,$5,$8}' \
       tmp/queryDB.lookup tmp/foldseek_output.txt > tmp/foldseek_log.txt

# Build tmalign
mkdir -p tmalign && cd tmalign 
curl https://zhanggroup.org/TM-align/TMalign.cpp > TMalign.cpp
patch -u TMalign.cpp -i ../tmalign_measurment.patch
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp 
cd ..

find pdbs/ -name 'd*' -printf "%f\n" > tmp/chain_list.txt

./tmalign/TMalign -dir ./pdbs/ tmp/chain_list.txt -fast \
  | tee tmp/tmalign_output.txt \
  | awk '/Name of Chain_1:/{printf "%s ",$4}
         /Name of Chain_2:/{printf "%s ",$4}
         /Length of Chain_1:/{printf "%s ",$4}
         /Length of Chain_2:/{printf "%s ",$4}
         /myduration/{printf "%s ", $2}
         /normalized by length of Chain_1/{printf "%s ",$2}
         /normalized by length of Chain_2/{printf "%s\n",$2}' \
  | awk '{sub(".*/", "", $3); sub(".*/", "", $2); print $2,$3,$4,$5,$6,$7,$1}' \
  > tmp/tmalign_log.txt

# Results
wc -l tmp/*_log.txt

