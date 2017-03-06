CC = gcc
CFLAGS += -fopenmp
LDFLAGS += -lm 

all: run

run: zeta2
	./zeta2
	
test: test0
	./test

zeta2: zeta2.c riemann.o
	${CC} zeta2.c ${CFLAGS} -o zeta2 riemann.o ${LDFLAGS}

test0: test.c riemann.o
	${CC} test.c ${CFLAGS} -o test riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c ${CFLAGS}
	
clean:
	-rm -f zeta2 test *.txt *.o
