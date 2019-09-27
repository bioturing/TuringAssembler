CC = gcc 

CXX = g++ 

CPP = cpp

LIBS = -pthread -O3 -std=c++11 \
       -Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
       -Llibs KMC/kmc_lib.a libs/libbz2.a libs/libz.a \
       libs/libbwa.a -lm 

# KMC_LIBS =  KMC/kmc_lib.a KMC/kmer_counter/libs/libz.a KMC/kmer_counter/libs/libbz2.a

GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -m64 -O3 -Wfatal-errors -Wall -Wextra \
         -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
         -DGIT_SHA='"$(GIT_SHA)"' \
         -Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
         -I ./src  -fPIC -g 

EXEC = skipping

EXEC_RELEASE = /src/skipping_static

# SRC = $(wildcard src/*.c)

SRC = src/assembly_graph.c 				\
      src/barcode_hash.c 				\
      src/barcode_builder.c 				\
      src/barcode_resolve2.c 				\
      src/basic_resolve.c 				\
      src/buffer_file_wrapper.c 	\
      src/dqueue.c 					\
      src/fastq_producer.c 				\
      src/get_buffer.c 					\
      src/io_utils.c 					\
      src/KMC_reader.c 					\
      src/kmer_build.c 					\
      src/kmhash.c 						\
      src/process.c 					\
      src/pthread_barrier.c 				\
      src/semaphore_wrapper.c 				\
      src/sort_read.c 						\
      src/time_utils.c 					\
      src/utils.c 					\
      src/verbose.c 					\
      src/build_bridge.c 				\
      src/read_list.c 					\
      src/map_contig.c 					\
      src/graph_search.c 				\
      src/helper.c 					\
      src/fastg.c 					\
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
      src/read_list.c 								\
      src/kmer_hash.c 								\
      src/fastq_reducer.c 						\
      src/main.c
      src/log.c

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

$(EXEC_RELEASE): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $(LDFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: debug
debug: CFLAGS += -fsanitize=address -fno-omit-frame-pointer -g  -gdwarf-3
debug: LIBS += -fsanitize=address -fno-omit-frame-pointer -lasan
debug: CC = gcc
debug: CXX = g++
debug: $(EXEC)


.PHONY: release
release: LIBS = -pthread -static -O3 -std=c++11 \
       -Wl,--whole-archive              \
       -lpthread KMC/libkmc.a \
       libs/libz.a libs/libbz2.a libs/libbwa.a \
       -Wl,--no-whole-archive -lm 
release: $(EXEC_RELEASE)

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ) $(EXEC) $(DEP)
