
# standard compile options for the c++ executable
GCC = gcc
EXE = simu
OBJS =  main.c
FLAGS = -O3
LIBFLAGS = -lm

# default super-target
all: $(EXE)

# the standard executable
$(EXE): $(OBJS)
	$(GCC) $(FLAGS) $^ -o $@ $(LIBFLAGS)

.PHONY: clean tidy

tidy:
	@find | egrep "#" | xargs rm -f
	@find | egrep "\~" | xargs rm -f
	@find | egrep ".txt" | xargs rm -f

clean: $(EXE)
	rm -f $(EXE)
