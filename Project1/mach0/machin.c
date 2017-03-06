#include <math.h>
#include "machin.h"

double arctan(double x, int n) {
	double sum = 0;
	
	for(int i = 1; i <= n; i++) {
		sum += pow(-1, i-1)*(pow(x, 2*i-1))/(2*i-1);
	}
	
	return sum;
}

double machin(int n) {
	return 4*(4*arctan(1.0/5, n)-arctan(1.0/239, n));
}