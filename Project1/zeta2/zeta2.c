#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "riemann.h"

int main ( int argc, char **argv ){
    
    int n;
	printf("Enter a value : ");
	scanf("%d", &n);
	
	double pi = riemann(n);
	
	printf("%f\n", pi);
}