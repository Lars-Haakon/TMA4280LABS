#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <mpi.h>

#include "machin.h"

double *global_vector;
double *local_vector;
double *partial_sums;

extern int my_rank, comm_sz;

double arctan(double x, int n, Method method) {
	if(my_rank == 0) {
		global_vector = malloc(n*sizeof(double));
		partial_sums = malloc(comm_sz*sizeof(double));

		// fill global vector
		for(int i = 1; i <= n; i++) {
			global_vector[i-1] = pow(-1, i-1)*(pow(x, 2*i-1))/(2*i-1);
		}
	}
	local_vector = malloc((n/comm_sz)*sizeof(double));

	MPI_Scatter(global_vector, n/comm_sz, MPI_DOUBLE, local_vector, n/comm_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double partial_sum = 0;
	for(int i = 0; i < n/comm_sz; i++) {
		partial_sum += local_vector[i];
	}
	free(local_vector);
	
	double sum = 0;
	if(method == ALLREDUCE) {
		MPI_Allreduce(&partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	} else if(method == RECURSIVEDOUBLE) {
		// Recursive doubling sum
		int power = 1;
		for(int d = 0; d < log2(comm_sz); d++) {
			
			int sigma = my_rank ^ power;

			double received_sum;
			MPI_Send(&partial_sum, 1, MPI_DOUBLE, sigma, 0, MPI_COMM_WORLD);
			MPI_Recv(&received_sum, 1, MPI_DOUBLE, sigma, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			partial_sum += received_sum;
			
			power = 2 << d;
		}
		sum = partial_sum;
	}
	
	if(my_rank == 0) {
		free(partial_sums);
		free(global_vector);
	}
	
	return sum;
}

double machin(int n, Method method) {
	
	return 4*(4*arctan(1.0/5, n, method) - arctan(1.0/239, n, method));
}