CC = gcc
LDFLAGS += -lm

all: run

run: mach0
	./mach0

utest: utest0
	./utest
	
vtest: vtest0
	./vtest
	
mach0: mach0.c machin.o
	${CC} mach0.c -o mach0 machin.o ${LDFLAGS}

utest0: utest.c machin.o
	${CC} utest.c -o utest machin.o ${LDFLAGS}
	
vtest0: vtest.c machin.o
	${CC} vtest.c -o vtest machin.o ${LDFLAGS}
	
machin.o: machin.c machin.h
	${CC} -c machin.c
	
clean:
	-rm -f mach0 utest vtest *.txt *.o
