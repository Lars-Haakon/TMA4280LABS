#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <unistd.h>

#include "riemann.h"

double walltime() {
    static struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec + 1e-6 * t.tv_usec);
}

int main ( int argc, char **argv ) {
	
	int p;
	printf("How many threads should be run? ");
	scanf("%d", &p);
	omp_set_num_threads(p);
	
	FILE* f = fopen("test.txt", "w");
	
	for(int k = 1; k <= 24; k++) {
		int n = 2 << (k-1);
		
		double start = walltime();
		double pi_computed = riemann(n);
		double finish = walltime();
		
		fprintf(f, "%d %e %.15f\n", k, finish-start, fabs(M_PI-pi_computed));
	
		sleep(1);
	}
	
	fclose(f);
}