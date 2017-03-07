#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "riemann.h"

double riemann(int n) {
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for(int i = 1; i <= n; i++) {
		//printf("Iteration: %d, Thread: %d of %d\n", i, omp_get_thread_num(), omp_get_num_threads());
		sum += 1/pow(i, 2);
	}
	
	return sqrt(6*sum);
}