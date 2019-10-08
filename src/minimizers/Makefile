CC = gcc

CXX = g++

CPP = cpp

LIBS = -pthread -std=c++11 \
       ../../libs/bzip2/libbz2.a ../../libs/zlib/libz.a \
       ../../libs/bwa/libbwa.a -lm


GIT_SHA := $(shell git rev-parse HEAD)

CFLAGS = -std=gnu99 -m64 -Wfatal-errors -Wall -Wextra \
         -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable \
         -DGIT_SHA='"$(GIT_SHA)"' \
         -Wl,--whole-archive -lpthread -Wl,--no-whole-archive \
         -I .. -g 

EXEC = skiminizer
SKIPPING_LIB = ../libskipping.a

SRC =  \
	   process.c \
	   opt.c \
	   count_barcodes.c \
	   smart_load.c \
	   ../dqueue.c \
	   ../fastq_producer.c \
	   ../utils.c \
	   ../io_utils.c \
	   ../get_buffer.c \
	   ../time_utils.c \
	   ../verbose.c \
	   ../pthread_barrier.c \
	   ../buffer_file_wrapper.c \
	   ../sort_read.c \
	   ../log.c

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
