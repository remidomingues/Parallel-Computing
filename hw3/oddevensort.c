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

double* merge_split(double* array, double* tmp, int length, int length2, int higher) {
    int i, start, end, step, idx, tmp_idx;
    double* result = (double*)malloc(length * sizeof(double));

    if(higher == 0) {
        start = 0;
        end = length;
        step = 1;
        idx = 0;
        tmp_idx = 0;
    } else {
        start = length-1;
        end = -1;
        step = -1;
        idx = length-1;
        tmp_idx = length2-1;
    }

    for(i = start; i != end; i += step) {
        if((array[idx] < tmp[tmp_idx] && higher == 0) || (array[idx] > tmp[tmp_idx] && higher == 1)) {
            result[i] = array[idx];
            idx += step;
        } else {
            result[i] = tmp[tmp_idx];
            tmp_idx += step;
        }
    }

    return result;
}

double* exchange(double* array, int length, int length2, int rank, int sendFirst) {
    double* tmp = (double*)malloc(length * sizeof(double));
    double* result;

    if(sendFirst == 1) {
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        MPI_Recv(tmp, length2, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result = merge_split(array, tmp, length, length2, 0);
    } else {
        MPI_Recv(tmp, length2, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(array, length, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        result = merge_split(array, tmp, length, length2, 1);
    }

    free(tmp);
    free(array);

    return result;
}

int getRemoteLength(rank, remoteRank, length, R) {
    if(R != 0 && (((R - rank <= 0) && (R - remoteRank > 0)) || ((R - rank > 0) && (R - remoteRank <= 0)))) {
        if(rank > remoteRank) {
            return length+1;
        }
        return length-1;
    }
    return length;
}

int main(int argc, char **argv)
{
    int i, j, step, I2;
    int P;
    int rank;
    clock_t begin, end;
    double time_spent;
    begin = clock();

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Find problem size N from command line */
    if (argc < 2) {
        if(rank == 0) {
            printf("No size N given\n");
        }
        return;
    }
    int N = atoi(argv[1]);
    /* local size */
    int R = N%P;
    int I = (N-R)/P; //Number of local elements common to all
    if((R-rank) > 0) {
        ++I;
    }

    double* x = (double*)malloc(I * sizeof(double));
    /* Initialization of the seed for the random generator */
    srandom(rank+1);

    /* data generation by each processor (avoid useless communication step for random data) */
    for (i = 0; i < I; i++) {
        x[i] = ((double)random())/RAND_MAX;
    }
    printf("%d > Local data at step %d is: ", rank, 0);
    display(x, I, rank);

    // Odd-even transposition sort
    for(step = 0; step < P; ++step) {
        // Efficient sequential sort in N log(N)
        qsort(x, I, sizeof(double), compare);
        // Even
        if(step % 2 == 0) {
            if(rank % 2 == 0) {
                if(rank < P-1) {
                    // Exchange with the next processor
                    I2 = getRemoteLength(rank, rank+1, I, R);
                    //printf(">>> %d > Step %d, echange with %d\n", rank, step, rank+1);
                    x = exchange(x, I, I2, rank+1, 1);
                }
            } else if(rank > 0) {
                // Exchange with the previous processor
                I2 = getRemoteLength(rank, rank-1, I, R);
                //printf(">>> %d > Step %d, echange with %d\n", rank, step, rank-1);
                x = exchange(x, I, I2, rank-1, 0);
            }
        // Odd
        } else {
            if(rank % 2 != 0) {
                if(rank < P-1) {
                    // Exchange with the next processor
                    I2 = getRemoteLength(rank, rank+1, I, R);
                    //printf(">>> %d > Step %d, echange with %d\n", rank, step, rank+1);
                    x = exchange(x, I, I2, rank+1, 1);
                }
            } else if(rank > 0) {
                // Exchange with the previous processor
                I2 = getRemoteLength(rank, rank-1, I, R);
                //printf(">>> %d > Step %d, echange with %d\n", rank, step, rank-1);
                x = exchange(x, I, I2, rank-1, 0);
            }
        }
    }

    // Data gathering
    double** data;
    if(rank == 0) {
        data = (double **)malloc(P * sizeof(double*));
    }

    if(rank == 0) {
        data[0] = x;
        for(i = 1; i < P; ++i) {
            if(R != 0 && R - i <= 0) {
                data[i] = (double *)malloc((I-1) * sizeof(double));
                MPI_Recv(data[i], I-1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                printf("%d > Data received from to %d: ", rank, i);
                display(data[i], I-1, rank);
            } else {
                data[i] = (double *)malloc(I * sizeof(double));
                MPI_Recv(data[i], I, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                printf("%d > Data received from to %d: ", rank, i);
                display(data[i], I, rank);
            }
        }
    } else {
        // printf("%d > Sending final data to %d: ", rank, 0);
        // display(x, I, rank);
        MPI_Send(x, I, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        free(x);
    }

    //Final display
    if(rank == 0) {
        printf("%d > Final result: ", 0);
        for(i = 0; i < P; ++i) {
            if(R != 0 && R - i <= 0) {
                step = I - 1;
            } else {
                step = I;
            }

            for(j = 0; j < step; j++) {
                printf("%.3f ", data[i][j]);
            }
        }
        printf("\n");
    }

    //Free memory
    if(rank == 0) {
        for(i = 0; i < P; ++i) {
            free(data[i]);
        }
        free(data);
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("0 > ELAPSED TIME: %.3f seconds", time_spent);
    }

    MPI_Finalize();
    return 0;
}
