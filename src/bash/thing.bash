#!/usr/bin/env bash

# uncomment to help with debugging
#set -o errexit
#set -o nounset
#set -o pipefail
#set -o xtrace

function test {

    cd ../../build
    make
    cd -
    cd ../../testsuite
    ./fast_test.sh
    cd -

    echo "applying age distributions..."

    mkdir -p out

    # parse the original ons household distribution into our own intermediate format
    parse-ons ../../data/ons-hh.csv out/tmp > out/distribution.csv
    
    # overlay the distribution on the output of ebola_build.exe, using an accept-reject algorithm
    accept-reject out/distribution.csv ../../testsuite/output/pop1_nodes.out out/tmp 10 > out/households.csv

    rm out/tmp* out/distribution.csv

    arcs-out-to-csv ../../testsuite/output/pop1_arcs.out > out/hh-connections.csv

    oversample out/households.csv > out/people.csv

    echo "=> out/households.csv, out/hh-connections.csv, out/people.csv"
}
function parse-ons {
    ons_csv="$1"
    shift

    sed "/,,,,/d" $ons_csv |
        sed "s/[(^M)|(a-z)|(A-Z)|(\ )]//g" |
        tail -n +2 |
        awk -F',' 'BEGIN{OFS=","}{if ($5!=0) print $5,$1,$2,$3,$4; }' |
        sort -nr |
        awk -F',' 'BEGIN{OFS=","; print "hh-size","a0-19","a20-64","a65+","freq" }{print $2,$3,$4,$5,$1;}'
}
function accept-reject {
    distribution="$1"
    shift

    nodes_out="$1"
    shift

    tmp_file="$1"
    shift

    max_hh_size="$1"
    shift

    echo "house-number","long","lat","a0-19","a20-64","a65+"

    sed "s/\t/,/g" $nodes_out > $tmp_file'_nodes'

    for i in `seq 1 $max_hh_size`
    do
        tail -n +2 $distribution | awk -F',' -v hhs=$i '{if ($1==hhs) print $n}' | sed "s/\r//g" > $tmp_file'_dist'
        awk -v hhs=$i -f accept-reject.awk $tmp_file'_dist' $tmp_file'_nodes'
    done
}
function oversample {
    households="$1"
    shift

    tail -n +2 $households | awk -L',' -f oversample.awk
}
function arcs-out-to-csv {
    nodes_out="$1"
    shift

    sed '/^$/d' ../../testsuite/output/pop1_arcs.out |
        awk 'BEGIN{OFS=","; print "hh-i,hh-j"}{ if(NR%2==1){p=$1}else{ print p,$1 }}'
}

# this calls the function provided as argument one,
# with subsequent arguments passed through
$@

