CC = mpicc
LDFLAGS += -lm

all: run

run: mach1
	mpirun -n 8 ./mach1
	
test: test0
	mpirun -n 8 ./test

mach1: mach1.c machin.o
	${CC} mach1.c -o mach1 machin.o ${LDFLAGS}

test0: test.c machin.o
	${CC} test.c -o test machin.o ${LDFLAGS}
	
machin.o: machin.c machin.h
	${CC} -c machin.c
	
clean:
	-rm -f mach1 test test.txt *.o
