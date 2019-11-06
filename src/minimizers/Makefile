CC = clang

CXX = clang++ 

CPP = cpp

LIBS = -pthread -O3 -std=c++11 \
       -Wl -lpthread \
       ../../libs/bzip2/libbz2.a ../../libs/zlib/libz.a \
       ../../libs/bwa/libbwa.a -lm 


GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -m64 -O3 -Wfatal-errors -Wall -Wextra \
         -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
         -DGIT_SHA='"$(GIT_SHA)"' \
         -Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
         -I .. -g 

EXEC = skiminizer

SRC =  \
	   process.c \
	   opt.c \
	   count_barcodes.c \
      ../buffer_file_wrapper.c 	\
      ../buffer_file_wrapper.c 	\
      ../dqueue.c 					\
      ../fastq_producer.c 				\
      ../log.c 						\
      ../get_buffer.c 					\
      ../io_utils.c 					\
      ../pthread_barrier.c 				\
      ../semaphore_wrapper.c 				\
      ../sort_read.c 						\
      ../time_utils.c 					\
      ../utils.c 					\
      ../verbose.c 					\

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

$(EXEC_RELEASE): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $(LDFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ) $(EXEC) $(DEP)
