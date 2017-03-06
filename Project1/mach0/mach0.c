#include <stdio.h>
#include <stdlib.h>

#include "machin.h"

int n;

int main ( int argc, char **argv ) {
	/*if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }
    n = atoi(argv[1]);*/
	int n;
	printf("Enter a value : ");
	scanf("%d", &n);
	
	double pi = machin(n);
	
	printf("%f\n", pi);
	
}