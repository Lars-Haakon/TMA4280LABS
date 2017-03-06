#include <math.h>

#include "riemann.h"

double riemann(int n) {
	double sum = 0;
	for(int i = 1; i <= n; i++) {
		sum += 1.0/(i*i);
	}
	
	return sqrt(6*sum);
}