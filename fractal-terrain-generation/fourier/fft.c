#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include <mpi.h>

#define TAU 6.28318530718 // see tauday.com

/* Check the return code of an MPI function */
void check(int rc);
/* Initialize the MPI communication system */
void initMPI(int argc, char **argv, int *N, int *rank);
/* Deals with command line arguments */
void parseCommandLine(int argc, char **argv, ptrdiff_t* N0, ptrdiff_t* N1, double* alpha) ; 
/* Generates the correct frequency space values */
void initialiseData(fftw_complex* data, ptrdiff_t N0, ptrdiff_t N1, double alpha, ptrdiff_t local_n0, int rank);
/* Outputs data to stdout */
void printData(fftw_complex* data, ptrdiff_t local_n0, ptrdiff_t N1, int rank, int N);

/* Computes a distributed fractal terrain generation using the Fourier tansform of pink noise 
 *
 * Compiling: mpicc SOURCE_PATH -lm -o BINARY_PATH
 * Execution: mpirun -np N ./BINARY_PATH width height p
 * With p: relaxation power, must be greater than 1. The higher it is, the smoother the terrain will be.
 */
int main(int argc, char** argv)
{
    int rank, N;
    clock_t begin, end;
    double alpha;

    ptrdiff_t N0, N1;
    parseCommandLine(argc, argv, &N0, &N1, &alpha);
    
    fftw_plan plan;
    fftw_complex *data;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;

    begin = clock();
    initMPI(argc, argv, &N, &rank);

    /* Initialise random generator seed */
    srandom(time(NULL) + rank);

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
            &local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);

    /* Create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
            FFTW_FORWARD, FFTW_ESTIMATE);

    /* Initialize data */
    initialiseData(data, N0, N1, alpha, local_n0, rank);

    /* compute transforms, in-place, as many times as desired */
    fftw_execute(plan);

    /* Calculate computation time */
    if(rank == 0) 
    {
        end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        fprintf(stderr,"0: Computation time %.3f seconds\n", time_spent);
    }

    /* Print results */
    printData(data, local_n0, N1, rank, N);

    /* Do FFTW and MPI cleanups */
    fftw_free(data);
    fftw_destroy_plan(plan);

    MPI_Finalize();

    return 0;
}

void parseCommandLine(int argc, char **argv, ptrdiff_t* N0, ptrdiff_t* N1, double* alpha) 
{
    int correct_argc = 3;
    if(argc < correct_argc)
    {
        fprintf(stderr,"Invalid number of arguments (%d), %d expected at least.\n", argc-1, correct_argc-1);
        fprintf(stderr,"Command line template is \"mpirun -np N ./BINARY_PATH width height p \n");
        exit(1);
    }

    *N0 = atoi(argv[1]);
    *N1 = atoi(argv[1]);
    *alpha = atof(argv[3]);
    if(*N0 <= 0 || *N1 <= 0)
    {
        fprintf(stderr,"Error: width and height should be positive\n");
        exit(1);
    }
    if(*alpha <= 0)
    {
        fprintf(stderr,"Error: for physical results, alpha should be positive.\n");
	exit(1);
    }
}

void check(int rc)
{
    if(rc != MPI_SUCCESS) 
    {
        char error_string[BUFSIZ];
        int error_class, string_length;
        MPI_Error_class(rc, &error_class);
        MPI_Error_string(rc, error_string, &string_length);
        fprintf(stderr,"MPI error %d: %s\n", error_class, error_string);
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(2);
    }
}

void initMPI(int argc, char **argv,int *N,int *rank) 
{
    check(MPI_Init(&argc, &argv));
    check(MPI_Comm_size(MPI_COMM_WORLD, N));
    check(MPI_Comm_rank(MPI_COMM_WORLD, rank));
}

void initialiseData(fftw_complex* data, ptrdiff_t N0, ptrdiff_t N1, double alpha, ptrdiff_t local_n0, int rank)
{
    int i,j;
    int fi;
    double r2,theta;
    for (i = 0; i < local_n0; ++i)
    {
        fi = local_n0*rank + i;
        if(fi < N0)
        {
            for (j = 0; j < N1/2; ++j)
            {
                /* Compute Box-Muller Transform */
                r2 = -2.0 * log((double)random()/RAND_MAX);
                theta = (double)random()/RAND_MAX;
                /* frequency space is now a normally distributed complex field */
                /* with a 1/f^alpha damping */
                data[i*N1+j] = pow(fi*fi + j*j,-alpha) * sqrt(r2) * cexp(I *TAU * theta);
            }
        }
        else
        {
            for (j = 0; j < N1/2; ++j)
            {
                data[i*N1+j] = 0;
            }
        }
        for (j = N1/2; j < N1; ++j)
        {
            data[i*N1+j] = 0;
        }
    }
    if(rank == 0) data[0] = 0;
}

void printData(fftw_complex* data, ptrdiff_t local_n0, ptrdiff_t N1, int rank, int N)
{
    int r,i,j;

    for(r = 0; r < N; r++)
    {
        if(r == rank)
        {
            for(i = 0; i < local_n0; i++)
            {
                printf("%f",creal(data[i*N1]));
                for(j = 1; j < N1; j++)
                {
                    printf(",%f",creal(data[i*N1+j]));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

