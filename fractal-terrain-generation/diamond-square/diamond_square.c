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

/*
 * Return an array containing the result of the diamond pass for the current iteration
 * on the right or down edge
 * data: the real sparse matrix
 * height: the virtual height
 * width: the virtual width
 * direction: RIGHT or DOWN
 */
float * getDataToSend(float ** data, int height, int width, int scalingFactor, direction_t direction) {
    float * result;
    int i, j, k = 0;

    if(direction == RIGHT) {
        result = (float *) malloc (height / 2 * sizeof(float));
        i = (height - 1) * scalingFactor;
        for(j = 1; j < width; j += 2, ++k) {
            result[k] = data[i][j * scalingFactor];
        }

    } else if(direction == DOWN) {
        result = (float *) malloc (width / 2 * sizeof(float));
        j = (width - 1) * scalingFactor;
        for(i = 1; i < height; i += 2, ++k) {
            result[k] = data[i * scalingFactor][j];
        }
    }
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

/* Diamond-square algorithm */
void diamondSquare(float ** data, int width, int height, int iterations) {
    int i, j, realI, realJ, k, scalingFactor, start, n, sum;

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

        // Diamond pass
        start = 1;
        // ==> Idée: on vire les IF de diamond pass (on calcule pas les bords)
        // Ensuite, à la fin de l'iteration, on échange les data, PUIS on calcule les diamants des bords
        // (attention si on est un process de bord)
        // et c'est là qu'on fait les if du coup
            // Du coup, on va devoir faire un red black est/west, north/south
            // les process enverront les C sur les bords et une fois reçus, les diamants des bords seront calculés
        for(i = 0; i < height; ++i) {
            // TODO: pas calculer les bords (changer idx start end)
            for(j = start; j < width; ++j) {
                n = 0;
                sum = 0;
                //TODO: keep for later
                if(i != 0) {
                    sum += data[realI-1][realJ];
                    ++n;
                }
                if(j != 0) {
                    sum += data[realI][realJ-1];
                    ++n;
                }
                if(i != 0) {
                    sum += data[realI+1][realJ];
                    ++n;
                }
                if(i != 0) {
                    sum += data[realI][realJ+1];
                    ++n;
                }
                //TODO: n = 4
                data[realI][realJ] = (sum + RANDOM_FACTOR * random() / (float) RAND_MAX) / n;
            }
            start = (start + 1) % 2;
        }
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
    diamondSquare(data, initWidth, initHeight, height, width, iterations);
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