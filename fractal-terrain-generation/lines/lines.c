#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

/* Check the return code of an MPI function */
void check(int rc);
/* Initialize the MPI communication system */
void initMPI(int argc, char **argv, int *N, int *rank);
/* Deals with command line arguments */
void parseCommandLine(int argc, char **argv, int* width, int* height, float* quality); 
/* calculate height for a certain number of iterations */
void calculateHeight(int width, int height, float quality, int N, int* data);
/* Outputs data to stdout */
void printData(int width, int height, int* data);
      
/* Computes a distributed fractal terrain using line displacement 
 *
 * Compiling: mpicc SOURCE_PATH -lm -o BINARY_PATH
 * Execution: mpirun -np N ./BINARY_PATH width height quality
 * quality: number of iterations is quality*width*height. Between 0 and 1.
 */

int main(int argc, char **argv)
{
    int  rank, N;
    clock_t begin, end;
    float time_spent;
    float quality;
    int width, height;

    parseCommandLine(argc, argv, &width, &height, &quality);

    begin = clock();
    initMPI(argc, argv, &N, &rank);

    /* Initialization of the random generator's seed */
    srandom(time(NULL) + rank);

    int* data = calloc(width*height,sizeof(int));
    assert(data != NULL);

    // Compute height matrix
    calculateHeight(width, height, quality, N,  data);
    MPI_Barrier(MPI_COMM_WORLD);
    check(MPI_Allreduce(MPI_IN_PLACE, data, width*height, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
    // Time display
    if(rank == 0) 
    {
        end = clock();
        time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
        fprintf(stderr,"0: Computation time %.3f seconds\n", time_spent);
        // Print results
        printData(width, height, data);
    }

    // Release memory
    free(data);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(int argc, char **argv, int* width, int* height, float* quality) 
{
    int correct_argc = 3;
    if(argc < correct_argc)
    {
        fprintf(stderr,"Invalid number of arguments (%d), %d expected at least.\n", argc-1, correct_argc-1);
        fprintf(stderr,"Command line template is \"mpirun -np N ./BINARY_PATH width height quality \n");
        exit(1);
    }

    *width = atoi(argv[1]);
    *height = atoi(argv[2]);
    *quality = atof(argv[3]);
    if(*width <= 0 || *height <= 0)
    {
        fprintf(stderr,"Error: width and height should be positive\n");
        exit(1);
    }
    if(*quality < 0 || *quality >1)
    {
        fprintf(stderr,"For physical results, %f should be between 0 and 1. The program may not terminate\n",*quality);
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

void calculateHeight(int width, int height, float quality, int N, int* data)
{
    const int Nlines = floor(quality*height*width/N);
    int i, x, y;
    float x1, y1, y2, x2;
    int sign;
    int ylim;

    for(i = 0; i < Nlines; i++)
    {
        x1 = (float)random()/RAND_MAX; y1 = (float)random()/RAND_MAX;
        x2 = (float)random()/RAND_MAX; y2 = (float)random()/RAND_MAX;
        sign = (random()%2)? -1 : 1;
        for(x = 0; x < width; x++)
        {
            ylim = floor((y2-y1)/(x1+x2)*(x - width*x1) + height*y1);
            if(ylim > height) ylim = height;
            for(y=0; y<ylim; y++)
            {
                data[height*x + y] += sign;
            }
        }
    }
}

void printData(int width, int height, int* data)
{
    int x,y;
    for(x = 0; x < width; x++)
    {
        printf("%d",data[height*x]);
        for(y = 1; y < height; y++)
        {
            printf(",%d",data[height*x + y]);
        }
        printf("\n");
    }
    fflush(stdout);
}
