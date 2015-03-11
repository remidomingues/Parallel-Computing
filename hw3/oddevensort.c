#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

void display(double* array, int length, int rank) {
    int i;
    for(i = 0; i < length; ++i) {
        printf("%d > %.2f ", rank, array[i]);
    }
    printf("\n");
}

int compare(const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}

double* merge_split(double* array, double* tmp, int length, int higher) {
    int i;
    int idx = 0;
    int tmp_idx = 0;
    double* result = (double*)malloc(length * sizeof(double));

    for(i = 0; i < length; ++i) {
        if((array[idx] < tmp[tmp_idx] && higher == 0) || (array[idx] > tmp[tmp_idx] && higher == 1)) {
            result[i] = array[idx];
            ++idx;
        } else {
            result[i] = tmp[tmp_idx];
            ++tmp_idx;
        }
    }
}

double* exchange(double* array, int length, int rank, int sendFirst) {
    double* tmp = (double*)malloc(length * sizeof(double));
    double* result;

    if(sendFirst == 1) {
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        MPI_Recv(tmp, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result = merge_split(array, tmp, length, 0);
    } else {
        MPI_Recv(tmp, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        result = merge_split(array, tmp, length, 1);
    }

    free(tmp);
    free(array);
}

int main(int argc, char **argv)
{
    int i, step;
    unsigned int P;
    unsigned int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Find problem size N from command line */
    if (argc < 2) {
        printf("No size N given\n");
        exit(-1);
    }
    int N = atoi(argv[1]);
    /* local size. Modify if P does not divide N */
    int I = N/P;
    double* x = (double*)malloc(I * sizeof(double));
    /* Initialization of the seed for the random generator */
    srandom(rank+1);

    /* data generation by each processor (avoid useless communication step for random data) */
    for (i = 0; i < I; i++) {
        x[i] = ((double)random())/RAND_MAX;
    }

    // Odd-even transposition sort
    unsigned int evenprocess = (i+1) % 2;
    unsigned int evenphase = 1;

    for(step = 0; step < N-1; ++step) {
        printf("%d > Internal sort:", rank);
        display(x, I*P, rank);

        // Efficient sequential sort in N log(N)
        qsort(x, I, sizeof(double), compare);

        printf("%d > become:", rank);
        display(x, I*P, rank);

        if((evenphase == 1 && evenprocess == 1) || (evenphase == 0 && evenprocess == 0)) {
            if(rank < N-1) {
                exchange(x, I, rank+1, 1);
            }
        } else {
            if(rank > 0) {
                exchange(x, I, rank-1, 0);
            }
        }

        evenphase = (evenphase + 1) % 2;
    }

    // Data gathering
    double* data;
    if(rank == 0) {
        data = (double *)malloc(P * I * sizeof(double));
    }

    MPI_Gather(x, I, MPI_DOUBLE, data, I, MPI_DOUBLE, 1, MPI_COMM_WORLD);

    //TODO REMOVE
    if(rank == 0) {
        display(data, I*P, rank);
    }
}
