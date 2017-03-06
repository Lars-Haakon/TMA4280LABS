CC = mpicc
CFLAGS += -fopenmp
LDFLAGS += -lm

all: run

run: zeta
	mpirun -n 8 ./zeta
	
test: test0
	mpirun -n 8 ./test

zeta: zeta.c riemann.o
	${CC} zeta.c ${CFLAGS} -o zeta riemann.o ${LDFLAGS}

test0: test.c riemann.o
	${CC} test.c ${CFLAGS} -o test riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c ${CFLAGS}
	
clean:
	-rm -f zeta test *.txt *.o
