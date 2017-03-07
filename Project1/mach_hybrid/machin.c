#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "machin.h"

double *global_vector;
double *local_vector;
double *partial_sums;

extern int my_rank, comm_sz;

double arctan(double x, int n) {
	if(my_rank == 0) {
		global_vector = malloc(n*sizeof(double));
		partial_sums = malloc(comm_sz*sizeof(double));

		// fill global vector
#pragma omp parallel for
		for(int i = 1; i <= n; i++) {
			//printf("Iteration: %d, Thread: %d of %d\n", i, omp_get_thread_num(), omp_get_num_threads());
			global_vector[i-1] = pow(-1, i-1)*(pow(x, 2*i-1))/(2*i-1);
		}
	}
	local_vector = malloc((n/comm_sz)*sizeof(double));

	MPI_Scatter(global_vector, n/comm_sz, MPI_DOUBLE, local_vector, n/comm_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double partial_sum = 0;
#pragma omp parallel for reduction(+:partial_sum)
	for(int i = 0; i < n/comm_sz; i++) {
		partial_sum += local_vector[i];
	}
	free(local_vector);
	
	MPI_Gather(&partial_sum, 1, MPI_DOUBLE, partial_sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(my_rank == 0) {
		double sum = 0;
#pragma omp parallel for reduction(+:sum)
		for(int i = 0; i < comm_sz; i++) {
			sum += partial_sums[i];
		}
		free(partial_sums);
		free(global_vector);
		
		return sum;
	}
	
	return 0;
}

double machin(int n) {
	
	return 4*(4*arctan(1.0/5, n) - arctan(1.0/239, n));
}