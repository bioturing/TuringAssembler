CC = gcc

CXX = g++

CPP = cpp

LIBS = -pthread -lm -lz -static -O3 -std=c++11 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

KMC_LIBS =  KMC/kmc_lib.a KMC/kmer_counter/libs/libz.a KMC/kmer_counter/libs/libbz2.a

GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -m64 -static -O3 -Wfatal-errors -Wall -Wextra \
	-Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
	-DGIT_SHA='"$(GIT_SHA)"' \
	-Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
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
      src/main.c

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ) $(KMC_LIBS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) $(KMC_LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $(LDFLAGS) $< -MM -MT $(@:.d=.o) >$@

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
