#include <stdio.h>
#include <mpi.h>

// Export a rendering matrix into an ASCII file.
// This one can be read in Matlab using "load -ascii rendering.txt"
void exportRendering(unsigned char * rendering, unsigned int height, unsigned int width) {
    FILE *fp = fopen("rendering.txt","w");
    int i, j;
    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++)
            fprintf(fp, "%hhu ", rendering[i+j*height]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

unsigned char compute_pixel(unsigned int d, unsigned int b, unsigned int N) {
    unsigned char value = 1;
    unsigned int z = 0;

    while((abs(z) < b) & (value < N)) {
        z = z^2+d;
        value = value+1;
    }

    return value;
}

void computeRendering(unsigned int P, unsigned int N, unsigned int b, unsigned int height,
        unsigned int width, unsigned char * rendering, unsigned int p) {
    unsigned int dx = 2*b/(width-1);
    unsigned int dy = 2*b/(height-1);
    unsigned int wp = width/P;
    unsigned int hp = height;
    unsigned int xoff = p*width/P;
    unsigned int yoff = 0;
    int dreal, dimag, d, x, y, i = 0;

    for(x = 0; x < wp-1; ++x) {
        dreal = (x+xoff)*dx-b;

        for(y = 0; y < height-1; ++y) {
            dimag = (y+yoff)*dy-b;
            d = dreal+i*dimag;

            rendering[i] = compute_pixel(d, b, N); //TODO CHECK i
            ++i; //TODO CHECK
        }
    }
}

int main(int argc, char **argv) {
    unsigned int N = 256;
    unsigned int b = 2;
    unsigned int P;
    unsigned int height = 2048;
    unsigned int width = 2048;
    unsigned int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned int size = height*width;
    unsigned int bufferSize = height*width/P;
    unsigned int dx = 2*b/(width-1);
    unsigned int dy = 2*b/(height-1);
    unsigned int wp = width/P;
    unsigned int hp = height;
    unsigned int xoff;
    unsigned int yoff = 0;
    int dreal, dimag, d, q, i, j;

    if(rank != 0) {
        size = bufferSize;
    }

    unsigned char rendering[size];
    computeRendering(P, N, b, height, width, rendering, rank);

    if(rank == 0) {
        MPI_Status status;
        unsigned char buffer[bufferSize];

        for(q = 1; q < P-1; ++q) {
            MPI_Recv(buffer, bufferSize, MPI_CHAR, q, 1, MPI_COMM_WORLD, &status);
            /*
            xoff = q*width/P;

            for(int x = 0; x < wp-1; ++x) {
                dreal = (x+xoff)*dx-b;

                for(int y = 0; y < height-1; ++y) {
                    dimag = (y+yoff)*dy-b;
                    d = dreal+i*dimag;

                    rendering[?] = buffer[?]; //TODO
                    ++i //TODO
                }
            }
            */

            for(i = q*wp, j = 0; i < (q+1)*wp-1; i += height, ++j) { //TODO: This is very likely to be wrong
                rendering[i] = buffer[j];
            }
        }
        exportRendering(rendering, height, width);
    } else {
        MPI_Send(rendering, size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
