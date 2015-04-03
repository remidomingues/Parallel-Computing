#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Direction */
typedef enum {UP, RIGHT, LEFT, DOWN} direction_t;

/* Default values */
WATER_HEIGHT = -1;
RANDOM_FACTOR = 10;
FACTOR_MULTIPLIER = 0.5;
MALTAB_OUTPUT = 1;

/* Initialize the values P and Q in order to distribute
 * N processes on a P*Q grid */
void getProcessesDistribution(int N, int *P, int *Q) {
    int i;

    *P = 1;
    *Q = N;

    for(i = 2; i < N; ++i) {
        if(i >= m2) {
            break;
        }

        if(N % i == 0) {
            *P = i;
            *Q = N / i;
        }
    }
}

/* Return the process rank given its coordinates on the grid */
int coordsToRank(int p, int q, int P) {
    return p * P + q;
}

/* Initialize the process' grids coordinates given its rank */
void rankToCoords(int r, int P, int * p, int * q) {
    *p = r / P;
    *q = r % P;
}

/* Return 1 if the grid coordinates are inside the grid boundaries
 * 0 else */
int existProcess(int p, int q, int P, int Q) {
    if(p >= 0 && p < P && q >= 0 && q < Q) {
        return 1;
    }
    return 0;
}

/* Return the final grid size (after n iterations) given n and the initial grid size */
int initialToFinalLength(int length, int iterations) {
    return (length - 1) * pow(2, iterations) - 1;
}

/* Return a scaling factor giving real i2 and j2 coordinates in
 * the data structure from virtual i and j coordinates
 * i2 = i * scalingFactor
 * j2 = j * scalingFactor */
int virtualToRealCoordsFactor(int iter, int iterations) {
    return pow(2, iterations - iter);
}

void parseCommandLine(int argc, char **argv, char ** inputFile, ** outputFile, int * iterations) {
    if (argc < 4) {
        printf("Invalid number of arguments (%d), %d expected.\n" +
            "Command line template is \"mpirun -np N ./BINARY_PATH FILE_PATH.in " +
            "FILE_PATH.out ITERATIONS [WATER_HEIGHT] [RAND_FACTOR] [FACTOR_MULTIPLIER]\"")
        exit(-1);
    }

    *inputFile = argv[1];
    *outputFile = argv[2];
    *iterations = atoi(argv[3])

    if(argc > 4) {
        WATER_HEIGHT = atoi(argv[4])
    }

    if(argc > 5) {
        RANDOM_FACTOR = atoi(argv[5])
    }

    if(argc > 6) {
        FACTOR_MULTIPLIER = atoi(argv[6])
    }
}

void initMPI(int argc, char **argv, int *N, int *rank) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, N);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

/* Return a sparse data matrix according to the input file and scaling factor */
float ** readData(char * inputFile, int * height, int * width, int * initWidth, int * initHeight, int iterations) {
    float ** data;
    int i, j, realI, realJ, scalingFactor;
    FILE * file = fopen(inputFile, "r");

    // Read initial data dimensions
    fscanf(file, "%d", initWidth);
    fscanf(file, "%d", initHeight);

    // Initialize final data matrix
    scalingFactor = virtualToRealCoordsFactor(0, iterations);
    *height = initialToFinalLength(*initHeight, iterations);
    *width = initialToFinalLength(*initWidth, iterations);
    data = (float **) malloc ((*height) * sizeof(float *));
    for(i = 0; i < height; ++i) {
        data[i] = (float *) malloc ((*width) * sizeof(float));
    }

    // Read file and fill the sparse final data matrix
    for(i = 0; i < initHeight; ++i) {
        realI = i * scalingFactor;
        for(j = 0; j < initWidth; ++j) {
            fscanf(file, "%f", &data[realI][j * scalingFactor]);
            printf("%f ", data[realI][j * scalingFactor]);
        }
        printf("\n")
    }
    fclose (file);

    return data;
}

int getDataLength(int rank, int length, int N) {
    int size = length / N;

    if(rank < length % N) {
        ++size;
    }

    if(rank != N - 1) {
        ++size;
    }

    return size;
}

float ** importAndScatterData(int rank, char * inputFile, int * height, int * width, int * initWidth, int * initHeight, int iterations, int P, int Q) {
    int r, i, j, k, realI, currentI = 0, currentJ = 0, tmpLength = -1, N = P*Q, p, q;
    int scalingFactor = virtualToRealCoordsFactor(0, iterations);
    float ** data, float * tmp;
    int * dataSize = (int *) malloc(2 * sizeof(int));

    if(rank == 0) {
        // Read data from file
        data = readData(inputFile, height, width, initWidth, initHeight, iterations);

        // Send data to each process
        for(r = 1; r < N; ++r) {
            rankToCoords(r, P, &p, &q);

            // Send data size
            dataSize[0] = getDataLength(r, initWidth, P);
            dataSize[1] = getDataLength(r, initHeight, Q);
            MPI_Send(dataSize, 2, MPI_INT, r, 1, MPI_COMM_WORLD);

            //Send data
            if(tmpLength != dataSize[0] * dataSize[1]) {
                tmpLength = dataSize[0] * dataSize[1];
                if(tmpLength != -1) {
                    free(tmp);
                }
                tmp = (float *) malloc (tmpLength * sizeof(float));
            }

            for(k = 0, i = currentI; i < currentI + dataSize[1]; ++i) {
                for(j = currentJ; j < currentJ + dataSize[0]; ++j, ++k) {
                    tmp[k] = data[i * scalingFactor][j * scalingFactor];
                }
            }

            MPI_Send(tmp, tmpLength, MPI_FLOAT, r, 1, MPI_COMM_WORLD);

            currentJ += dataSize[0];
            currentI += dataSize[1];
        }

    } else {
        // Receive data size
        MPI_Recv(dataSize, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *initHeight = dataSize[0];
        *initWidth = dataSize[1];
        *height = initialToFinalLength(*initHeight, iterations);
        *width = initialToFinalLength(*initWidth, iterations);

        // Receive data
        tmpLength = (*initHeight) * (*initWidth);
        tmp = (float *) malloc (tmpLength * sizeof(float));
        MPI_Recv(tmp, tmpLength, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Initialize data matrix
        data = (float **) malloc ((*height) * sizeof(float *));
        for(i = 0; i < height; ++i) {
            data[i] = (float *) malloc ((*width) * sizeof(float));
        }

        // Fill data matrix
        for(k = 0, i = 0; i < *initHeight; ++i) {
            realI = i * scalingFactor;
            for(j = 0; j < *initWidth; ++j, ++k) {
                data[realI][j * scalingFactor] = tmp[k];
            }
        }
    }

    // Release memory
    if(tmpLength != -1) {
        free(tmp);
    }
    free(dataSize);
}

/*
 * Return an array containing the result of the diamond pass for the current iteration
 * on the right or down edge
 * data: the real sparse matrix
 * height: the virtual height
 * width: the virtual width
 * If the direction is RIGHT or DOWN, the random factors are computed and added to the squares to send
 */
float * getDataToSend(float ** data, int height, int width, int scalingFactor, direction_t direction) {
    float * result;
    int i, j, k = 0;

    if(direction == TOP) {
        result = (float *) malloc ((width - 1) / 2 * sizeof(float));
        i = scalingFactor;
        for(j = 1; j < width; j += 2, ++k) {
            result[k] = data[i][j * scalingFactor];
        }

    } else if(direction == RIGHT) {
        result = (float *) malloc ((height - 1) / 2 * sizeof(float));
        j = (width - 2) * scalingFactor;
        for(i = 1; i < height; i += 2, ++k) {
            result[k] = data[i * scalingFactor][j] + RANDOM_FACTOR * random() / (float) RAND_MAX;
        }

    } else if(direction == DOWN) {
        result = (float *) malloc ((width - 1) / 2 * sizeof(float));
        i = (height - 2) * scalingFactor;
        for(j = 1; j < width; j += 2, ++k) {
            result[k] = data[i][j * scalingFactor] + RANDOM_FACTOR * random() / (float) RAND_MAX;
        }

    } else if(direction == LEFT) {
        result = (float *) malloc ((height - 1) / 2 * sizeof(float));
        j = scalingFactor;
        for(i = 1; i < height; i += 2, ++k) {
            result[k] = data[i * scalingFactor][j];
        }
    }

    return result;
}

void updateDiamonds(float ** data, int height, int width, int scalingFactor, direction_t direction, float * squares, float * remoteSquares) {
    int i, realI, realJ, k = 0, n;

    if(direction == RIGHT || direction == LEFT) {
        if(direction == RIGHT) {
            realJ = (width - 1) * scalingFactor;
        } else if(direction == LEFT) {
            realJ = scalingFactor;
        }

        for(i = 1; i < height; i += 2, ++k) {
            n = 2;
            realI = i * scalingFactor;
            data[realI][realJ] = data[realI-1][realJ] + squares[k] + data[realI+1][realJ];
            if(remoteSquares != NULL) {
                data[realI][realJ] = (data[realI][realJ] + remoteSquares[k]) / 4.0
            } else {
                data[realI][realJ] /= 3.0;
            }
        }

    } else if(direction == TOP || direction == DOWN) {
        if(direction == TOP) {
            realI = scalingFactor;
        } else if(direction == DOWN) {
            realI = (height - 1) * scalingFactor;
        }

        for(j = 1; j < width; j += 2, ++k) {
            realJ = J * scalingFactor;
            data[realI][realJ] = data[realI][realJ-1] + squares[k] + data[realI][realJ+1];
            if(remoteSquares != NULL) {
                data[realI][realJ] = (data[realI][realJ] + remoteSquares[k]) / 4.0
            } else {
                data[realI][realJ] /= 3.0;
            }
        }
    }
}

/* Exchange squares and update diamonds for one specific direction
 * p and q are the coordinates of the REMOTE process */
void exchangeSquares(float ** data, int height, int width, int scalingFactor, int p, int q, int P, int Q,
    direction_t direction, int sendFirst) {
    float * squares, * remoteSquares;
    int length = height;

    if(direction == TOP || direction == DOWN) {
        length = width;
    }

    squares = getDataToSend(data, height, width, scalingFactor, direction);
    if(existProcess(p, q, P, Q) {
        if(sendFirst == 1) {
            MPI_Send(squares, (length - 1) / 2, MPI_FLOAT, coordsToRank(p, q, P), 1, MPI_COMM_WORLD);
            remoteSquares = (float *) malloc ((length - 1) / 2 * sizeof(float));
            MPI_Recv(remoteSquares, (length - 1) / 2, MPI_FLOAT, coordsToRank(p, q, P), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            remoteSquares = (float *) malloc ((length - 1) / 2 * sizeof(float));
            MPI_Recv(remoteSquares, (length - 1) / 2, MPI_FLOAT, coordsToRank(p, q, P), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(squares, (length - 1) / 2, MPI_FLOAT, coordsToRank(p, q, P), 1, MPI_COMM_WORLD);
        }

        updateDiamonds(data, height, width, scalingFactor, direction, squares, remoteSquares);
        free(remoteSquares);

    } else {
        updateDiamonds(data, height, width, scalingFactor, direction, squares, NULL);
    }
    free(squares);
}

/* Send and receive the squares on the four sides according to the red black algorithm
 * Diamond pass is completed after this exchange */
void exchangeSquares(float ** data, int height, int width, int scalingFactor, int p, int q, int P, int Q) {
    // Right - left exchange
    if(p % 2 == 0) {
        exchangeSquares(data, height, width, scalingFactor, p+1, q, P, Q, RIGHT, 1);
        exchangeSquares(data, height, width, scalingFactor, p-1, q, P, Q, LEFT, 1);
    } else {
        exchangeSquares(data, height, width, scalingFactor, p-1, q, P, Q, LEFT, 0);
        exchangeSquares(data, height, width, scalingFactor, p+1, q, P, Q, RIGHT, 0);
    }

    // Top - down exchange
    if(q % 2 == 0) {
        exchangeSquares(data, height, width, scalingFactor, p, q+1, P, Q, DOWN, 1);
        exchangeSquares(data, height, width, scalingFactor, p, q-1, P, Q, TOP, 1);
    } else {
        exchangeSquares(data, height, width, scalingFactor, p, q-1, P, Q, TOP, 0);
        exchangeSquares(data, height, width, scalingFactor, p, q+1, P, Q, DOWN, 0);
    }
}

/* Diamond-square algorithm */
void diamondSquare(float ** data, int width, int height, int iterations, int p, int q, int P, int Q) {
    int i, j, realI, realJ, k, scalingFactor, start;

    // Iterations
    for(k = 0; k < iterations; ++k) {
        scalingFactor = virtualToRealCoordsFactor(k, iterations);
        width = 2 * width - 1;
        height = 2 * height - 1;

        // Square pass
        for(i = 1; i < height; i += 2) {
            realI = i * scalingFactor;
            for(j = 1; i < width; j += 2) {
                realJ = j * scalingFactor;
                data[realI][realJ] = (data[realI-1][realJ-1] + data[realI-1][realJ+1] +
                    data[realI+1][realJ-1] + data[realI+1][realJ+1] + RANDOM_FACTOR * random() / (float) RAND_MAX) / 4.0;
            }
        }

        // Diamond pass - Diamonds on the side are not computed
        start = 1;
        for(i = 1; i < height - 1; ++i) {
            realI = i * scalingFactor;
            for(j = start + 1; j < width - 1; ++j) {
                realJ = j * scalingFactor;
                data[realI][realJ] = (data[realI-1][realJ] + data[realI][realJ-1] +
                    data[realI+1][realJ] + data[realI][realJ+1] + RANDOM_FACTOR * random() / (float) RAND_MAX) / 4.0;
            }
            start = (start + 1) % 2;
        }

        // Data exchange - Only the squares on the side and the random factors are exchanged
        exchangeSquares(data, height, width, scalingFactor, p, q, P, Q);

        // Random factor attenuation
        RANDOM_FACTOR *= FACTOR_MULTIPLIER;
    }
}

void pourWater(float ** data, int height, int width) {
    int i, j;
    for(i = 0; i < height; ++i) {
        for(j = 0; j < width; ++j) {
            if(data[i][j] < WATER_HEIGHT) {
                data[i][j] = WATER_HEIGHT;
            }
        }
    }
}

/* Compute a distributed fractal terrain generation using the diamond-square algorithm
 * Command line is:
 * mpirun -np N ./BINARY_PATH FILE_PATH.in FILE_PATH.out ITERATIONS [WATER_HEIGHT] [RAND_FACTOR] [FACTOR_MULTIPLIER]
 * With N: number of processes, FILE_PATH.in the input file, FILE_PATH.out the output file
 * ITERATIONS: the scaling factor between the input and output data structure (> 0)
 *
 * To display the content of the output file, use
 * surf(0:100:(WIDTH-1)*100,0:100:(HEIGHT-1)*100, TERRAIN);
 */
int main(int argc, char **argv)
{
    int rank, N, p, P, q, Q;
    int heigth, width, initHeight, initWidth, iterations;
    clock_t begin, end;
    double time_spent;
    char *inputFile, *outputFile;
    float **data;

    parseCommandLine(argc, argv, &inputFile, &outputFile, &iterations);
    initMPI(argc, argv, &N, &rank)
    getProcessesDistribution(N, &P, &Q);
    rankToCoords(rank, P, &p, &q)

    /* Initialization of the random generator's seed */
    srandom(time(NULL) + rank);

    begin = clock();
    data = importAndScatterData(rank, inputFile, &height, &width, &initWidth, &initHeight, iterations, P, Q);
    diamondSquare(data, initWidth, initHeight, height, width, iterations, p, q, P, Q);
    pourWater(data, height, width);
    //TODO gather
    //TODO export

    // Release memory
    //TODO

    //TODO TEST FUNCTIONS

    // Time display
    if(rank == 0) {
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("0 > Computation time: %.3f seconds", time_spent);
    }

    MPI_Finalize();
    return 0;
}
