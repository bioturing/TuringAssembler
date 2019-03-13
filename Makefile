CC = gcc

CPP = cpp

LIBS = -pthread -flto -lm -lz

CFLAGS = -Wfatal-errors -Wall -Wextra -fPIC -std=gnu99 -g

EXEC = skipping

SRC = $(wildcard src/*.c)

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

-include $(dep)

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

