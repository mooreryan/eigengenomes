#!/bin/bash

make clean && make

kmer_len=3
num_hps=4

rm test_files/s?.hash.*

echo;echo "Running hash_and_count"
parallel "echo;echo;valgrind ./hash_and_count 0 $kmer_len $num_hps {} > {.}.hash.counts" ::: test_files/*fastq

diff <(sort test_files/s1.hash.counts) <(sort hash_and_count.s1.expected)
if [ $? -ne 0 ]
then
    echo "s1 outfile (left) doesn't match expected (right)!"
    exit 1
fi

diff <(sort test_files/s2.hash.counts) <(sort hash_and_count.s2.expected)
if [ $? -ne 0 ]
then
    echo "s2 outfile (left) doesn't match expected (right)!"
    exit 1
fi

echo;echo "Running weight_counts"
valgrind ./weight_counts $num_hps test_files/*.hash.counts > test_out

diff test_out test_weight_counts.expected_output

if [ $? -ne 0 ]
then
    echo
    echo "The output hash changed!"
    diff -y test_out test_weight_counts.expected_output
    exit 1
else
    echo
    echo "It is all good!"
fi

rm test_out
