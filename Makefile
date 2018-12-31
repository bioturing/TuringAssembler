GCC = gcc

AR = ar

LIBS = -pthread -flto -lm -lz
#-fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wall -Wextra -fPIC -std=gnu99 -O2 -g

EXEC = skipping

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
      src/time_utils.c 				\
      src/utils.c 				\
      src/verbose.c 				\
      src/main.c

$(EXEC):
	$(CC) -o $(EXEC) $(CFLAGS) $(SRC) $(LIBS)

all: $(EXEC)

lazy: CFLAGS += -DUSE_PRIME_HASH
lazy: $(EXEC)

clean:
	rm -f $(EXEC)
