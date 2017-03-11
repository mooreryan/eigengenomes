#!/bin/bash

make clean && make

kmer_len=11
num_hps=11

dir=$PWD/test_files/from_lsa
input_files=$dir/*.fastq

echo "Input files: $input_files"

echo;echo "Running hash_and_count"
parallel "./hash_and_count 0 $kmer_len $num_hps {} > {}.hash_counts" \
         ::: $dir/*.fastq

if [ $? -ne 0 ]
then
    echo "Failed!" && exit
fi

echo;echo "Running weight_counts"
./weight_counts $num_hps $dir/*.fastq.hash_counts | \
    ./map_mm > thingy.mm 2> /dev/null

if [ $? -ne 0 ]
then
    echo "Failed!" && exit
fi


echo;echo "Running find_cluster_seeds.py"
LARGE=f python find_cluster_seeds.py $dir thingy.mm

if [ $? -ne 0 ]
then
    echo "Failed!" && exit
fi

rm thingy* $dir/*.hash_counts*
