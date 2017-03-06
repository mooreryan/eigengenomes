CC = gcc
CFLAGS = -Wall -g -O2
LSA_FLAGS = -lz -lm
OBJS = eigg_encode.o eigg_hash.o eigg_hyperplane.o eigg_kmers.o \
       eigg_print.o vendor/tommyarray.o vendor/tommyhashlin.o

.PHONY: all
all: hash_and_count weight_counts

hash_and_count: $(OBJS)
	$(CC) $(CFLAGS) $(LSA_FLAGS) -o $@ $^ $@.c

weight_counts:
	$(CC) $(CFLAGS) -lm -o $@ $@.c

vendor/tommyarray:
	$(CC) $(CFLAGS) -c -o $@.o $@.c

vendor/tommyhashlin:
	$(CC) $(CFLAGS) -c -o $@.o $@.c

tommy_test: vendor/tommyarray.o vendor/tommyhashlin.o
	$(CC) $(CFLAGS) -o $@ $^ $@.c

.PHONY: clean
clean:
	-rm hash_and_count $(OBJS)
	-rm -r hash_and_count.dSYM
	-rm weight_counts
	-rm -r weight_counts.dSYM
