#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "machin.h"

double arctan(double x, int n) {
	double sum = 0;

#pragma omp parallel for reduction(+:sum)
	for(int i = 1; i <= n; i++) {
		//printf("Iteration: %d, Thread: %d of %d\n", i, omp_get_thread_num(), omp_get_num_threads());
		sum += pow(-1, i-1)*(pow(x, 2*i-1))/(2*i-1);
	}
	
	return sum;
}

double machin(int n) {
	return 4*(4*arctan(1.0/5, n)-arctan(1.0/239, n));
}