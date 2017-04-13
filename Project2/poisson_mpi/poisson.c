/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;
int rank;
int size;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real** b2, size_t m, size_t n);
real rhs(real x, real y);
real exact_solution(real x, real y);
void export_plot_data(int m, real *grid, real **b);
void print_matrix(real** b, size_t m, size_t n);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    int n = atoi(argv[1]);
    int m = n/size;
    real h = 1.0 / n;
	
	if((n & (n - 1)) != 0 || (size & (size - 1)) != 0) {
		if(rank==0) {
			printf("both problem size and number of MPI procs has to be a power of 2\n");
		}
		exit(-1);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double start;
	if(rank==0) {
		printf("number of processes: %d\n", size);
		start = MPI_Wtime();
	}

	MPI_Datatype col, coltype;
	MPI_Type_vector(m, 1, n, MPI_DOUBLE, &col);       
	MPI_Type_commit(&col);
	MPI_Type_create_resized(col, 0, sizeof(double), &coltype);
	MPI_Type_commit(&coltype);
	
	int scount[size];
	int sdisp[size];
	int rcount[size];
	int rdisp[size];
	for(int i = 0; i < size; i++) {
		scount[i]=m;
		sdisp[i]=i*m;
		rcount[i]=m*m;
		rdisp[i]=i*m*m;
	}
	
    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(n, false);
    for (size_t i = 0; i < n; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    /*
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    real **b = mk_2D_array(m, n, false);
	real **b2 = mk_2D_array(m, n, false);
    real **bt = mk_2D_array(m, n, false);

    /*
     * This vector will holds coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doings 
     * reallocations at each function call.
     */
    int nn = 4 * n;
    real *z = mk_1D_array(nn, false);

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     * 
     */
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            b[i][j] = h * h * rhs(grid[rank*m+i+1], grid[j+1]);
        }
    }

    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for 
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input 
     * array (first argument) so that the initial values are overwritten.
     */
	 
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z, &nn);
    }
	
	MPI_Alltoallv(&(b[0][0]), scount, sdisp, coltype, &(b2[0][0]), rcount, rdisp, MPI_DOUBLE, MPI_COMM_WORLD);
	transpose(bt, b2, m, n);
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }
	
    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            bt[i][j] /= (diag[i+rank*m] + diag[j]);
        }
    }
	
    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
	MPI_Alltoallv(&(bt[0][0]), scount, sdisp, coltype, &(b2[0][0]), rcount, rdisp, MPI_DOUBLE, MPI_COMM_WORLD);
	transpose(b, b2, m, n);
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) {
		double finish = MPI_Wtime();
		
		printf("Elapsed time: %e seconds\n", finish-start);
	}
	
    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
    real e_local = 0.0;
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n-1; j++) {
			real exact = exact_solution(grid[rank*m+i+1], grid[j+1]);
			
            e_local = e_local > fabs(b[i][j]-exact) ? e_local : fabs(b[i][j]-exact);
        }
    }
	
	real e;
	MPI_Reduce(&e_local, &e, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(rank==0) {
		printf("e_max = %e\n", e);
	}
	MPI_Finalize();

    return 0;
}

void print_matrix(real** b, size_t m, size_t n) {
	if(rank!=0) {
		int r;
		MPI_Recv(&r, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			printf("%.15f ", b[i][j]);
		}
		printf("\n");
	}
	
	if(rank != size-1) {
		sleep(2);
		int s = 0;
		MPI_Send(&s, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) {
    //return 2 * (y - y*y + x - x*x);
	return 5*PI*PI*sin(PI*x)*sin(2*PI*y);
	//return exp(x)*sin(2*PI*x)*sin(2*PI*y);
	//return 1;
}

real exact_solution(real x, real y) {
	return sin(PI*x)*sin(2*PI*y);
}

void export_plot_data(int m, real *grid, real **b) {
	FILE* f = fopen("plotdata.txt", "w");

	fprintf(f, "%d\n", m);
	for(int i = 0; i < m; i++) {
		fprintf(f, "%.15f\n", grid[i+1]);
	}
	
	for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
			
			fprintf(f, "%.15f ", b[i][j]);
		}
		fprintf(f, "\n");
	}
	
	fclose(f);
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void transpose(real **bt, real** b2, size_t m, size_t n) {
   int n_ass = 0;
   for(int i = 0; i < m; i++) {
	   for(int j = 0; j < n; j+=m) {
		   real* t1 = &b2[i][j];
		   real* t2 = &bt[n_ass%m][m*(n_ass/m)];
		   
		   for(int k = 0; k < m; k++) {
			   t2[k] = t1[k];
		   }
		   n_ass++;
	   }
   }
}

/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}
