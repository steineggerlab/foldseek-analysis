#!/bin/bash
GAPOPEN=$1
GAPEXTEND=$2

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    exit
fi

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}

task(){
    tmp/ssw_test -o "$GAPOPEN" -e "$GAPEXTEND" -a tmp/s.mat -p \
        -f 50 tmp/target.fasta "$1" \
        2> /dev/null \
    | mawk '/^target/{target=$2}
            /^query/{query=$2}
            /^optimal_alignment_score/{score=$2; print query,target,score}' \
    | sort -k1,1 -k3,3nr \
    >  "tmp/alignments/${1##*/}"".m8";
}


cp tmp/sub_score.mat tmp/s.mat  # must be shorter than 16 characters

N=64
open_sem $N
for thing in tmp/splits/*; do
    run_with_lock task "$thing"
done

(while [[ $(pidof ssw_test) ]]; do sleep 1; done)

