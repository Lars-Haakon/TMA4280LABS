#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "riemann.h"

int main ( int argc, char **argv ) {
	double pi_expected = 3.1415926535897932;
	
	FILE* f = fopen("test.txt", "w");
	
	for(int k = 1; k <= 24; k++) {
		int n = 2 << (k-1);
		
		double pi_computed = riemann(n);
		
		fprintf(f, "%d %.20f %.20f\n", k, pi_computed, fabs(pi_expected-pi_computed));
	}
	
	fclose(f);
}