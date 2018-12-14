GCC = gcc

AR = ar

LIBS = -pthread -flto -lm -lz -fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wall -Wextra -fPIC -std=gnu99 -O2 -g

EXEC = kmer_count

# SRC = src/dqueue.c 				\
#       src/get_buffer.c 				\
#       src/graph.c 				\
#       src/io_utils.c 				\
#       src/kmer_count.c 				\
#       src/kmhash.c 				\
#       src/semaphore_wrapper.c 			\
#       src/verbose.c 				\
#       src/utils.c 				\
#       src/main.c

SRC = src/dqueue.c 				\
      src/fastq_producer.c 			\
      src/get_buffer.c 				\
      src/graph_assembly.c 			\
      src/io_utils.c 				\
      src/kmer_count.c 				\
      src/kmhash.c 				\
      src/semaphore_wrapper.c 			\
      src/verbose.c 				\
      src/utils.c 				\
      src/main.c

all:
	$(CC) -o $(EXEC) $(CFLAGS) $(SRC) $(LIBS)

clean:
	rm -f $(EXEC)
