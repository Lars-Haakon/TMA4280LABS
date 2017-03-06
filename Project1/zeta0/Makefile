CC = gcc
LDFLAGS += -lm

all: run

run: zeta0
	./zeta0

utest: utest0
	./utest
	
vtest: vtest0
	./vtest
	
zeta0: zeta0.c riemann.o
	${CC} zeta0.c -o zeta0 riemann.o ${LDFLAGS}

utest0: utest.c riemann.o
	${CC} utest.c -o utest riemann.o ${LDFLAGS}
	
vtest0: vtest.c riemann.o
	${CC} vtest.c -o vtest riemann.o ${LDFLAGS}
	
riemann.o: riemann.c riemann.h
	${CC} -c riemann.c
	
clean:
	-rm -f zeta0 utest vtest verification.txt *.o
