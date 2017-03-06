CC = mpicc
LDFLAGS += -lm

all: run

run: zeta1
	mpirun -n 2 ./zeta1
	
test: test0
	mpirun -n 2 ./test

zeta1: zeta1.c riemann.o
	${CC} zeta1.c -o zeta1 riemann.o ${LDFLAGS}

test0: test.c riemann.o
	${CC} test.c -o test riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c
	
clean:
	-rm -f zeta1 test test.txt *.o
