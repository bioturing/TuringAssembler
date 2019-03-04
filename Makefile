GCC = gcc

AR = ar

LIBS = -L./ -pthread -flto -lm -lz -fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wall -Wextra -fPIC -std=gnu99 -g

EXEC = skipping

SRC = 						\
      src/assembly_graph.c 			\
      src/dqueue.c 				\
      src/fastq_producer.c 			\
      src/get_buffer.c 				\
      src/io_utils.c 				\
      src/k63_build.c 				\
      src/k63_count.c 				\
      src/k63hash.c 				\
      src/k31_build.c 				\
      src/k31_count.c 				\
      src/k31hash.c 				\
      src/pthread_barrier.c 			\
      src/semaphore_wrapper.c 			\
      src/time_utils.c 				\
      src/utils.c 				\
      src/verbose.c 				\
      src/test_hash_count.c 			\
      src/main.c

OBJS= src/dqueue.o src/assembly_graph.o src/fastq_producer.o src/get_buffer.o \
	  src/io_utils.o src/k63_build.o src/k63_count.o src/k63hash.o src/k31_build.o \
	  src/k31_count.o src/k31hash.o src/pthread_barrier.o src/semaphore_wrapper.o \
	  src/time_utils.o src/utils.o src/verbose.o  src/test_hash_count.o

$(EXEC):
	$(CC) -o $(EXEC) $(CFLAGS) $(SRC) $(LIBS)

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

libskip.a: $(OBJS)
	$(AR) -csru $@ $(OBJS)

simple: $(OBJS) libskip.a src/simple.o
	$(CC) -o simple src/simple.o libskip.a $(CFLAGS) $(LIBS) 

all: $(EXEC)

time: $(EXEC)

mem: CFLAGS += -DUSE_BINARY_SEARCH
mem: $(EXEC)

clean:
	rm -f $(EXEC)
	rm -rf libskip.a
	rm -rf src/*.o
