CC = gcc

CPP = cpp

LIBS = -pthread -flto -lm -lz 

GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -O2 -Wfatal-errors -Wall -Wextra \
	-Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
	-DGIT_SHA='"$(GIT_SHA)"' \
	-I ./src \
	-g

EXEC = skipping


# SRC = $(wildcard src/*.c)

SRC = src/assembly_graph.c 				\
      src/barcode_retriever.c 				\
      src/barcode_resolve.c 				\
      src/basic_resolve.c 				\
      src/dqueue.c 					\
      src/fastq_producer.c 				\
      src/get_buffer.c 					\
      src/io_utils.c 					\
      src/k31_build.c 					\
      src/k31hash.c 					\
      src/k63_build.c 					\
      src/k63hash.c 					\
      src/kmer_count.c 					\
      src/process.c 					\
      src/pthread_barrier.c 				\
      src/semaphore_wrapper.c 				\
      src/test_hash_count.c 				\
      src/time_utils.c 					\
      src/utils.c 					\
      src/verbose.c 					\
      src/scaffolding/bin_hash.c 						\
      src/scaffolding/algorithm.c 						\
      src/scaffolding/compare.c 						\
      src/scaffolding/buck.c 						\
      src/scaffolding/contig.c 						\
      src/scaffolding/edge.c 						\
      src/scaffolding/global_params.c 						\
      src/scaffolding/output.c 						\
      src/scaffolding/score.c 						\
      src/scaffolding/scaffolding.c 						\
      src/scaffolding/scaffold.c 						\
      src/scaffolding/contig_graph.c 						\
      src/main.c

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CC) $(INCLUDE) $(LDFLAGS) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $(INCLUDE) $(LDFLAGS) $< -MM -MT $(@:.d=.o) >$@

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
