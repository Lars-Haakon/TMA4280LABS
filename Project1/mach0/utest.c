#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "machin.h"

int main ( int argc, char **argv ) {
	
	int n = 3;
	
	double pi_expected = 3.1415926535897932;
	double pi_computed = machin(n);
	
	printf("PI expected: %f\nPI computed: %f\nDifference: %f\n", 
			pi_expected, pi_computed, fabs(pi_expected-pi_computed));
}