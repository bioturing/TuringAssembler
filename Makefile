CC = gcc 

CXX = g++ 

CPP = cpp

LIBS = -pthread -std=c++11 \
       -Wl,--whole-archive -lpthread -Wl,--no-whole-archive libs/KMC/libkmc.a libs/zlib/libz.a libs/bzip2/libbz2.a  \
       libs/bwa/libbwa.a -lm

GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -m64 -O3 -Wfatal-errors -Wall -Wextra \
         -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
         -DLOG_USE_COLOR -DGIT_SHA='"$(GIT_SHA)"' \
         -Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
         -I ./ -I ./src  -fPIC -g -I ./src/minimizers \
	 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -fPIC

EXEC = TuringAssembler_dev

EXEC_RELEASE = TuringAssembler

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
      src/complex_resolve.c 				\
      src/split_molecules.c 				\
      src/cluster_molecules.c 				\
      src/simple_queue.c 				\
      src/khash_operations.c 				\
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
      src/read_list.c 								\
      src/kmer_hash.c 								\
      src/fastq_reducer.c 						\
      src/main.c \
      src/log.c 					\
	  src/build_hash_table.c\
	  src/resolve_big.c \
      src/unit_test.c \
      src/minimizers/minimizers.c \
      src/minimizers/smart_load.c \
      src/minimizers/count_barcodes.c

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
       -lpthread libs/KMC/libkmc.a \
       -Wl,--no-whole-archive -lm \
       libs/zlib/libz.a libs/bwa/libbwa.a \
       -Wl,--no-whole-archive -lm libs/bzip2/libbz2.a
release: CC = gcc 
release: CXX = g++
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
