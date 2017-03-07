#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <mpi.h>

#include "riemann.h"

double *global_vector;
double *local_vector;
double *partial_sums;

extern int my_rank, comm_sz;

double riemann(int n) {
	
	if(my_rank == 0) {
		global_vector = malloc(n*sizeof(double));
		partial_sums = malloc(comm_sz*sizeof(double));
		
		// fill global vector
		for(int i = 1; i <= n; i++) {
			global_vector[i-1] = 1/pow(i, 2);
		}
	}
	local_vector = malloc((n/comm_sz)*sizeof(double));
	
	MPI_Scatter(global_vector, n/comm_sz, MPI_DOUBLE, local_vector, n/comm_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double partial_sum = 0;
	for(int i = 0; i < n/comm_sz; i++) {
		partial_sum += local_vector[i];
	}
	free(local_vector);
	
	MPI_Gather(&partial_sum, 1, MPI_DOUBLE, partial_sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(my_rank == 0) {
		double sum = 0;
		for(int i = 0; i < comm_sz; i++) {
			sum += partial_sums[i];
		}
		
		free(partial_sums);
		free(global_vector);
		
		return sqrt(6*sum);
	}
	
	return 0;
}