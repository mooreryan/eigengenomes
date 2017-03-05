# Eigengenomes

A C implementation of Brian Cleary's [LSA](https://github.com/brian-cleary/LatentStrainAnalysis), for eigengenome partitioning.

## Install

Run these commands.

```bash
git clone https://github.com/mooreryan/eigengenomes.git
cd eigengenomes
make
```

This will make the executable file `eigg`. You can move this to somewhere on your path now. E.g.,

```bash
sudo mv eigg /usr/local/bin
```

### Uninstall

In the source directory, run `make clean`.

## Usage

The `eigg` program takes four command line arguments, and writes to standard out.

```
Usage: ./eigg <1: seed> <2: kmer size> <3: number of hyperplanes> <4: seqs.fastq> > out.hashq
```

## Error codes

- 0: Success
- 1: Argument error
- 2: Couldn't open a file
- 3: Not an even number of forward and reverse reads
- 4: A sequnece is shorter than the kmer length
