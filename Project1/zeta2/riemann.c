#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "riemann.h"

double riemann(int n) {
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for(int i = 1; i <= n; i++) {
		sum += 1/pow(i, 2);
	}
	
	return sqrt(6*sum);
}