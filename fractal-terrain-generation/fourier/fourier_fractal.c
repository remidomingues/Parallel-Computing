#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "mpi.h"

#define TAU 6.28318530718 // see tauday.com
#define RATIO 0.25 // ratio between frequency size and real size

/* Global width & height values */
int width, height;

/* Check the return code of an MPI function */
void check(int rc);
/* Initialize the MPI communication system */
void initMPI(int argc, char **argv, int *N, int *rank);
/* Deals with command line arguments */
void parseCommandLine(int argc, char **argv, double* p); 
/* Generates and shares a full array of complex coefficients */
void generateP(int rank, int N, double alpha, double complex* P);
/* Calculates the height with an inverse Fourier Transform */
void calculateHeight(double complex* P, int rank, int N, double* data);
/* Outputs data to stdout */
void printData(double* data, int rank, int N);
      
/* Computes a distributed fractal terrain generation using the Fourier tansform of pink noise 
 *
 * Compiling: mpicc SOURCE_PATH -lm -o BINARY_PATH
 * Execution: mpirun -np N ./BINARY_PATH width height p
 * With p: relaxation power, must be greater than 1. The higher it is, the smoother the terrain will be.
 */

int main(int argc, char **argv)
{
    int  rank, N;
    clock_t begin, end;
    double time_spent;
    double alpha;

    parseCommandLine(argc, argv, &alpha);

    begin = clock();
    initMPI(argc, argv, &N, &rank);

    /* Initialization of the random generator's seed */
    srandom(time(NULL) + rank);

    double complex* P = calloc(RATIO*RATIO*width*height,sizeof(double complex));
    assert(P != NULL);
    double* data = calloc(width*height/N,sizeof(double));
    assert(data != NULL);

    // Get full array of coefficients for each dimension
    generateP(rank, N, alpha, P);

    // Compute height sub-matrix
    calculateHeight(P, rank, N, data);

    // Time display
    if(rank == 0) 
    {
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        fprintf(stderr,"0: Computation time %.3f seconds\n", time_spent);
    }

    // Print results
    printData(data, rank, N);

    // Release memory
    free(P);
    free(data);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(int argc, char **argv, double* alpha) 
{
    int correct_argc = 3;
    if(argc < correct_argc)
    {
        fprintf(stderr,"Invalid number of arguments (%d), %d expected at least.\n", argc-1, correct_argc-1);
        fprintf(stderr,"Command line template is \"mpirun -np N ./BINARY_PATH width height p \n");
        exit(1);
    }

    width = atoi(argv[1]);
    height = atoi(argv[2]);
    *alpha = atof(argv[3]);
    if(width <= 0 || height <= 0)
    {
        fprintf(stderr,"Error: width and height should be positive\n");
        exit(1);
    }
    if(*alpha < 1)
    {
        fprintf(stderr,"For physical results, %f should be positive. Continuing...\n",*alpha);
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

void generateP(int rank, int N, double alpha, double complex* P)
{
    int x,y,fx;
    int Nx = RATIO*height;
    int K = RATIO*width/N;
    double damp = -0.5 * alpha;
    double r2,theta;
    double complex* Ptemp = calloc(K*Nx,sizeof(double complex));
    for(x = 0; x < K; x++)
    {
        fx = rank*K + x;
        for(y = 0; y < Nx; y++)
        {
            r2 = -2 * log(random()/(double)RAND_MAX);
            theta = random()/(double)RAND_MAX;
            Ptemp[Nx*x + y] = pow(fx*fx + y*y,damp) * sqrt(r2) * cexp(I * TAU * theta);

        }
    }
    if(rank == 0) { Ptemp[0] = 0; }
    check(MPI_Allgather(Ptemp,K*Nx,MPI_DOUBLE_COMPLEX,P,K*Nx,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD));
    free(Ptemp);

}

void calculateHeight(double complex* P, int rank, int N, double* data)
{
    int i,j,x,fx,fy;
    int Nx = RATIO*width;
    int Ny = RATIO*height;
    int K = width/N;
    double theta;

    for(i = 0; i < K; i++)
    {
        x = rank*K + i;
        for(j = 0; j < height; j++)
        {
            for(fx = 0; fx < Nx; fx++)
            {
                for(fy = 0; fy < Ny; fy++)
                {
                    theta = fx * x/(double)width + fy * j/(double)height;
                    data[height*i + j] +=  creal(P[fx*Nx + fy] * cexp(I * TAU * theta));
                }
            }
        }
    }
}

void printData(double* data, int rank, int N)
{
    int r,i,j;
    int K = width/N;

    for(r = 0; r < N; r++)
    {
        if(r == rank)
        {
            for(i = 0; i < K; i++)
            {
                printf("%f",data[height*i]);
                for(j = 1; j < height; j++)
                {
                    printf(",%f",data[height*i + j]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
