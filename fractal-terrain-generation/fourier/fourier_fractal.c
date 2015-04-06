#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "mpi.h"

#define TAU 6.28318530718 // see tauday.com

/* Global width & height values */
int width, height;

/* Check the return code of an MPI function */
void check(int rc);
/* Initialize the MPI communication system */
void initMPI(int argc, char **argv, int *N, int *rank);
/* Deals with command line arguments */
void parseCommandLine(int argc, char **argv, double* p); 
/* Generates and shares a full array of complex coefficients */
void generateP(int rank, double k, double p, double complex* P);
/* Calculates the height with an inverse Fourier Transform */
void calculateHeight(double complex* Px, double complex* Py, double alpha, int rank, int N, double* data);
/* Outputs data to stdout */
void printData(double* data, int rank, int N);
      
/* Computes a distributed fractal terrain generation using the Fourier tansform of pink noise 
 *
 * Compiling: mpicc diamond_square.c -lm -o diamond_square
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

    double complex* Px = calloc(width/2, sizeof(double complex));
    assert(Px != NULL);
    double complex* Py = calloc(height/2, sizeof(double complex));
    assert(Py != NULL);
    double* data = calloc(width*height/N,sizeof(double));
    assert(data != NULL);

    // Get full array of coefficients for each dimension
    generateP(rank, 0.5*width/N, alpha, Px);
    generateP(rank, 0.5*height/N, alpha, Py);

    calculateHeight(Px, Py, alpha, rank, N, data);

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
    free(Px);
    free(Py);
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

void generateP(int rank, double k, double alpha, double complex* P)
{
    int f;
    double r2,theta;
    double complex* Ptemp = calloc(k,sizeof(double complex));
    for(f = 0; f < k; f++)
    {
        r2 = -2 * log((double)rand()/(double)RAND_MAX);
        theta = (double)rand()/(double)RAND_MAX;
        Ptemp[f] = pow(rank*k+f+1,-alpha) * sqrt(r2) * cexp(I * TAU * theta);
    }
    check(MPI_Allgather(Ptemp,k,MPI_DOUBLE_COMPLEX,P,k,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD));
    free(Ptemp);
}

void calculateHeight(double complex* Px, double complex* Py, double alpha, int rank, int N, double* data)
{
    int i,j,fx,fy;
    int K = width/N;
    double complex theta;

    for(i = 0; i < K; i++)
    {
        for(j = 0; j < height; j++)
        {
            for(fx = 0; fx < 0.5*width; fx++)
            {
                for(fy = 0; fy < 0.5*height; fy++)
                {
//                     damping = pow(pow(fx+1.0, 2.0) + pow(fy+1.0, 2.0), -alpha/2.0);
                    theta = (fx+1.0) * (rank*K + i)/(double)width + (fy+1.0) * j/(double)height;
                    data[height*i + j] +=  creal(Px[fx] * Py[fy] * cexp(I * TAU * theta));
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
