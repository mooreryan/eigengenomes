#!/bin/bash

make clean && make

kmer_len=3
num_hps=4

rm test_files/s?.hash.*

echo;echo "Running hash_and_count"
parallel "echo;echo;valgrind ./hash_and_count 0 $kmer_len $num_hps {} > {.}.hash.counts" ::: test_files/s?.fastq

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
valgrind ./weight_counts $num_hps test_files/s?.hash.counts > test_out.mm

diff test_out.mm test_weight_counts.expected_output

if [ $? -ne 0 ]
then
    echo
    echo "The output hash changed!"
    diff test_out.mm test_weight_counts.expected_output
    exit 1
fi

echo;echo "Running map_mm"
valgrind ./map_mm < test_out.mm 1> test_out.mapped.mm 2> test_out.hb_map

diff test_out.mapped.mm test_out.mapped.mm.expected
if [ $? -ne 0 ]
then
    echo
    echo "The output of map_mm hash changed!"
    exit 1
fi

# grep out the valgrind output
diff <(grep -v "==" test_out.hb_map) test_out.hb_map.expected
if [ $? -ne 0 ]
then
    echo
    echo "The stderr of map_mm hash changed!"
    exit 1
fi


echo;echo "Running lsi.py"
python lsi.py test_files test_out.mapped.mm test_files/s?.hash.counts
if [ $? -ne 0 ]
then
    echo "Something went wrong in the lsi.py step"
    exit 1
fi

# diff test_files/kmer_lsi.gensim test_files/kmer_lsi.gensim.expected
# if [ $? -ne 0 ]
# then
#     echo "The kmer_lsi.gensim file doesn't match!"
#     exit 1
# fi

# diff test_files/kmer_lsi.gensim.projection test_files/kmer_lsi.gensim.projection.expected
# if [ $? -ne 0 ]
# then
#     echo "The kmer_lsi.gensim.projection file doesn't match!"
#     exit 1
# fi

echo;echo "It's all good!"
echo;echo

rm test_out.mm test_out.mapped.mm test_out.hb_map test_files/kmer_lsi.gensim.projection test_files/kmer_lsi.gensim
