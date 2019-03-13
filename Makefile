CC = gcc

CPP = cpp

LIBS = -pthread -flto -lm -lz

CFLAGS = -Wfatal-errors -Wall -Wextra -Wno-unused-function -fPIC -std=gnu99 -g

EXEC = skipping

SRC = $(wildcard src/*.c)

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ) $(EXEC) $(DEP)
