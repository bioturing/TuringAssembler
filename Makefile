GCC = gcc

AR = ar

LIBS = -pthread -flto -lm -lz -fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wall -Wextra -fPIC -std=gnu99 -O2 -g

EXEC = skipping

SRC = src/assembly_graph.c 			\
      src/dqueue.c 				\
      src/fastq_producer.c 			\
      src/get_buffer.c 				\
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

time: $(EXEC)

mem: CFLAGS += -DUSE_BINARY_SEARCH
mem: $(EXEC)

clean:
	rm -f $(EXEC)
