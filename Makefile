CC = gcc
CFLAGS = -Wall -g -O2
LSA_FLAGS = -lz -lm
OBJS = eigg_encode.o eigg_hash.o eigg_hyperplane.o eigg_kmers.o \
       eigg_print.o

.PHONY: all
all: hash_and_count weight_counts

hash_and_count: $(OBJS)
	$(CC) $(CFLAGS) $(LSA_FLAGS) -o $@ $^ $@.c

weight_counts:
	$(CC) $(CFLAGS) -lm -o $@ $@.c

.PHONY: clean
clean:
	-rm hash_and_count $(OBJS)
	-rm -r hash_and_count.dSYM
	-rm weight_counts
	-rm -r weight_counts.dSYM
