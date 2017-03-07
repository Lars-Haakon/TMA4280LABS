#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "machin.h"

double walltime() {
    static struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec + 1e-6 * t.tv_usec);
}

int main ( int argc, char **argv ) {
	FILE* f = fopen("test.txt", "w");
	
	for(int k = 1; k <= 24; k++) {
		int n = 2 << (k-1);
		
		double start = walltime();
		double pi_computed = machin(n);
		double finish = walltime();
		
		fprintf(f, "Elapsed time: %e seconds\n%d %.20f %.20f\n", finish-start, k, pi_computed, fabs(M_PI-pi_computed));
	}
	
	fclose(f);
}