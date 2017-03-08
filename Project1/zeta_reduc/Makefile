CC = mpicc
LDFLAGS += -lm

all: run

run: zeta
	mpirun -n 2 ./zeta
	
test1: test-a
	mpirun -n 8 ./test1
	
test2: test-r
	mpirun -n 8 ./test2

zeta: zeta.c riemann.o
	${CC} zeta.c -o zeta riemann.o ${LDFLAGS}

test-a: test1.c riemann.o
	${CC} test1.c -o test1 riemann.o ${LDFLAGS}

test-r: test2.c riemann.o
	${CC} test2.c -o test2 riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c
	
clean:
	-rm -f zeta test1 test2 *.txt *.o
