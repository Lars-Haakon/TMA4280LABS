#include <math.h>

#include "riemann.h"

double riemann(int n) {
	double sum = 0;
	for(int i = 1; i <= n; i++) {
		sum += 1/pow(i, 2);
	}
	
	return sqrt(6*sum);
}