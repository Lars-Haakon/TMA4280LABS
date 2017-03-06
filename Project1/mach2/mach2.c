#include <stdio.h>
#include <stdlib.h>

#include "machin.h"

int main ( int argc, char **argv ) {
	
	int n;
	printf("Enter a value : ");
	scanf("%d", &n);
	
	double pi = machin(n);
	
	printf("%f\n", pi);
	
}