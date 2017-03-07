#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>

#include "riemann.h"

int my_rank, comm_sz;

int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	assert(comm_sz%2==0);
	
	int n;
	if(my_rank == 0) {
		printf("Enter a value :\n");
		scanf("%d", &n);
		for(int i = 1; i < comm_sz; i++)
			MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	
	double pi_computed = riemann(n);

	if(my_rank == 0) {
		double finish = MPI_Wtime();
		
		printf("pi_error: %.20f\nElapsed time: %e seconds\n", fabs(M_PI-pi_computed), finish-start);
	}
	
	MPI_Finalize();
}