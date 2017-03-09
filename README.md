# Eigengenomes

A C and Python implementation of Brian Cleary's [LSA](https://github.com/brian-cleary/LatentStrainAnalysis), for eigengenome partitioning.

## Install

Run these commands.

```bash
git clone https://github.com/mooreryan/eigengenomes.git
cd eigengenomes
make
```

This will make the executable files needed for the pipeline.

### Uninstall

In the source directory, run `make clean`. Also, if you moved any of the program binaries, you will need to manually remove them.

## Usage

Run `hash_and_count` for each fastq file. It assumes reads are paired.

```
Usage: ./hash_and_count <1: seed> <2: kmer size> <3: number of hyperplanes> <4: seqs.fastq> > seqs.hash_counts
```

Run `weight_counts` once on all the output files from the `hash_and_count` program.

```
Usage: ./weight_counts <1: num hyperplanes used in previous step> *.hash_bucket_counts > cool_sample_group.hash_counts
```

## Example

There aren't wrapper scripts yet, but here is an example of running the steps implemented thus far.

`test_files` contains `s1.fastq` and `s2.fastq`.

This command will hash the kmers and count them. It will create these files: `test_files/s{1,2}.hash_counts`.

```bash
parallel "./hash_and_count 0 3 4 {} > {.}.hash_counts" ::: test_files/*fastq
```

This command will weight the counts using something similar to `tf-idf` weighting. It will create the file `test_files/all.hash_counts`.

```bash
./weight_counts 4 test_files/*.hash_counts > test_files/all.hash_counts
```

Aaaand, the remaining steps aren't finished yet ;)

## Error codes

- 0: Success
- 1: Argument error
- 2: Couldn't open a file
- 3: Not an even number of forward and reverse reads
- 4: A sequnece is shorter than the kmer length
- 5: A hashed kmer bucket was missing from a counting hash

## Issues

### TODO

- The counting arrays are way to memory intensive, consider hash tables for counting.
- Nucleotides are being mapped to 5 different integers instead of complex numbers.

### Differences from the original progam

- Hyperplanes are drawn at random through the origin rather than based on existing kmers.
- The idf score is calculated a bit differently.
- It is not designed to be run across many nodes of a compute cluster. This may change.
