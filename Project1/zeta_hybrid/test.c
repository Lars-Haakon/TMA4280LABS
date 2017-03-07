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
	
	FILE* f = fopen("test.txt", "w");
	
	for(int k = 3; k <= 24; k++) {
		int n = 2 << (k-1);
		
		MPI_Barrier(MPI_COMM_WORLD);
		double start = MPI_Wtime();
		
		double pi_computed = riemann(n);
		
		if(my_rank == 0) {
			double finish = MPI_Wtime();
			
			fprintf(f, "Elapsed time: %e seconds\n%d %.20f %.20f\n", finish-start, k, pi_computed, fabs(M_PI-pi_computed));
		}
	}
	
	fclose(f);
	
	MPI_Finalize();
}