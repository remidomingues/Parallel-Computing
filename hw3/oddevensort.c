#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

void display(double* array, int length, int rank) {
    int i;
    for(i = 0; i < length; ++i) {
        printf("%.2f ", array[i]);
    }
    printf("\n");
}

int compare(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
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

double* exchange(double* array, int length, int length2, int rank, int sendFirst) {
    double* tmp = (double*)malloc(length * sizeof(double));
    double* result;

    if(sendFirst == 1) {
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        MPI_Recv(tmp, length2, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result = merge_split(array, tmp, length, 0);
    } else {
        MPI_Recv(tmp, length2, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        result = merge_split(array, tmp, length, 1);
    }

    free(tmp);
    free(array);
}

int main(int argc, char **argv)
{
    int i, j, step, I2;
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
    /* local size */
    int R = N%P;
    int I = (N-R)/P; //Number of local elements common to all
    if(R - rank > 0) {
        ++I;
    }

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

    for(step = 0; step < P; ++step) {
        printf("%d > Internal sort:", rank);
        display(x, I*P, rank);

        // Efficient sequential sort in N log(N)
        qsort(x, I, sizeof(double), compare);

        printf("%d > become:", rank);
        display(x, I*P, rank);

        if((evenphase == 1 && evenprocess == 1) || (evenphase == 0 && evenprocess == 0)) {
            if(rank < N-1) {
                I2 = I;
                if(R != 0 && ((R - rank <= 0) != (R - rank + 1 <= 0))) {
                    I2 = I - 1;
                }
                exchange(x, I, I2, rank+1, 1);
            }
        } else {
            if(rank > 0) {
                I2 = I;
                if(R != 0 && ((R - rank <= 0) != (R - rank - 1 <= 0))) {
                    I2 = I + 1;
                }
                exchange(x, I, I2, rank-1, 0);
            }
        }

        evenphase = (evenphase + 1) % 2;
    }

    // Data gathering
    double** data;
    if(rank == 0) {
        data = (double **)malloc(P * sizeof(double*));
    }

    if(rank == 0) {
        for(i = 1; i < P; ++i) {
            if(R != 0 && R - i <= 0) {
                data[i] = (double *)malloc((I-1) * sizeof(double));
                MPI_Recv(data[i], I-1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                data[i] = (double *)malloc(I * sizeof(double));
                MPI_Recv(data[i], I, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    } else {
        MPI_Send(x, I, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    //Final display
    if(rank == 0) {
        printf("%d > Final result\n", 0);
        for(i = 0; i < P; ++i) {
            if(R != 0 && R - i <= 0) {
                step = I - 1;
            } else {
                step = I;
            }

            for(j = 0; j < step; j++) {
                printf("%.2f ", data[i][j]);
            }
        }
        printf("\n");
    }
}
