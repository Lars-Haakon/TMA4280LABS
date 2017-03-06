CC = gcc
CFLAGS += -fopenmp
LDFLAGS += -lm

all: run

run: mach2
	./mach2

test: test0
	./test
	
mach2: mach2.c machin.o
	${CC} mach2.c ${CFLAGS} -o mach2 machin.o ${LDFLAGS}

test0: test.c machin.o
	${CC} test.c ${CFLAGS} -o test machin.o ${LDFLAGS}
	
machin.o: machin.c machin.h
	${CC} -c machin.c ${CFLAGS}
	
clean:
	-rm -f mach2 test *.txt *.o
