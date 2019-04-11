CC = gcc

CPP = cpp

LIBS = -L/mnt/hdd1/tools/longline/htslib -pthread -flto -lm -lz -lhts

CFLAGS = -Wfatal-errors -Wall -Wextra -Wno-unused-function \
	-Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -fPIC -std=gnu99 -g

EXEC = skipping

# SRC = $(wildcard src/*.c)

SRC = src/assembly_graph.c 				\
      src/barcode_retriever.c 				\
      src/basic_resolve.c 				\
      src/dqueue.c 					\
      src/fastq_producer.c 				\
      src/get_buffer.c 					\
      src/io_utils.c 					\
      src/k31_build.c 					\
      src/k31_count.c 					\
      src/k31hash.c 					\
      src/k63_build.c 					\
      src/k63_count.c 					\
      src/k63hash.c 					\
      src/process.c 					\
      src/pthread_barrier.c 				\
      src/semaphore_wrapper.c 				\
      src/test_hash_count.c 				\
      src/time_utils.c 					\
      src/utils.c 					\
      src/verbose.c 					\
      src/main.c

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: debug
debug: LIBS += -fsanitize=undefined,address
debug: $(EXEC)

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ) $(EXEC) $(DEP)
