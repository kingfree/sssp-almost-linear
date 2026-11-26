CC = gcc
CFLAGS = -O2 -std=c11 -Wall -Wextra -Wno-unused-parameter
OBJS = util.o graph.o route.o owner.o outer.o hier.o theta.o eval.o main.o
sssp: $(OBJS)
	$(CC) $(CFLAGS) -o sssp $(OBJS) -lm
%.o: %.c sssp.h len.h
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f *.o sssp
