CC = mpicc
LDFLAGS += -lm

all: run

run: zeta
	mpirun -n 2 ./zeta
	
test: test0
	mpirun -n 8 ./test

zeta: zeta.c riemann.o
	${CC} zeta.c -o zeta riemann.o ${LDFLAGS}

test0: test.c riemann.o
	${CC} test.c -o test riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c
	
clean:
	-rm -f zeta test *.txt *.o
