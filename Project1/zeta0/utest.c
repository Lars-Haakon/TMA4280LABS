#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "riemann.h"

int main ( int argc, char **argv ) {
	
	int n = 3;
	
	double pi_computed = riemann(n);
	
	printf("PI expected: %f\nPI computed: %f\nDifference: %f\n", 
			M_PI, pi_computed, fabs(M_PI-pi_computed));
}