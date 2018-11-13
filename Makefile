GCC = gcc

AR = ar

LIBS = -pthread -flto -lm -lz

CFLAGS = -Wfatal-errors -Wall -fPIC -std=gnu99 -O2 -g

EXEC = kmer_count

SRC = src/dqueue.c 				\
      src/get_buffer.c 				\
      src/graph.c 				\
      src/io_utils.c 				\
      src/kmer_count.c 				\
      src/kmhash.c 				\
      src/verbose.c 				\
      src/utils.c 				\
      src/main.c

all:
	$(CC) -o $(EXEC) $(CFLAGS) $(SRC) $(LIBS)

clean:
	rm -f $(EXEC)
