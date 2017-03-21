CC = mpicc
LDFLAGS += -lm

all: run

run: mach
	mpirun -n 2 ./mach
	
test1: test-a
	mpirun -n 8 ./test1
	
test2: test-r
	mpirun -n 8 ./test2

mach: mach.c machin.o
	${CC} mach.c -o mach machin.o ${LDFLAGS}

test-a: test1.c machin.o
	${CC} test1.c -o test1 machin.o ${LDFLAGS}

test-r: test2.c machin.o
	${CC} test2.c -o test2 machin.o ${LDFLAGS}
	
machin.o: machin.c machin.h
	${CC} -c machin.c
	
clean:
	-rm -f mach test1 test2 *.txt *.o
