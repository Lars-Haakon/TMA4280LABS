CC = mpicc
LDFLAGS += -lm

all: run

run: mach
	mpirun -n 8 ./mach
	
test: test0
	mpirun -n 8 ./test

mach: mach.c machin.o
	${CC} mach.c -o mach machin.o ${LDFLAGS}

test0: test.c machin.o
	${CC} test.c -o test machin.o ${LDFLAGS}
	
machin.o: machin.c machin.h
	${CC} -c machin.c
	
clean:
	-rm -f mach test test.txt *.o
