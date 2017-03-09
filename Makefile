CC = gcc
CFLAGS = -Wall -g -O2
LSA_FLAGS = -lz -lm
OBJS = eigg_encode.o eigg_hash.o eigg_hyperplane.o eigg_kmers.o \
       eigg_print.o
TOMMY_OBJS = vendor/tommyarray.o vendor/tommyhashlin.o vendor/tommylist.o

.PHONY: all
all: hash_and_count weight_counts map_mm

hash_and_count: $(OBJS) $(TOMMY_OBJS)
	$(CC) $(CFLAGS) $(LSA_FLAGS) -o $@ $^ $@.c

make_mm:
	$(CC) $(CFLAGS) -o $@ $@.c

weight_counts: $(TOMMY_OBJS)
	$(CC) $(CFLAGS) -lm -o $@ $^ $@.c

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
