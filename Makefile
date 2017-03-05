CC = gcc
CFLAGS = -Wall -g -O2
LSA_FLAGS = -lz -lm
OBJS = eigg_encode.o eigg_hash.o eigg_hyperplane.o eigg_kmers.o \
       eigg_print.o

eigg : $(OBJS)
	$(CC) $(CFLAGS) $(LSA_FLAGS) -o eigg $(OBJS) eigg.c

.PHONY : clean
clean :
	-rm eigg eigg.dSYM $(OBJS)
	-rm -r eigg.dSYM
