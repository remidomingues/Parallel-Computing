#include <stdio.h>
#include <complex.h>
#include <mpi.h>

// Export a rendering matrix into an ASCII file.
// This one can be read in Matlab using "load -ascii rendering.txt"
void exportRendering(unsigned char * rendering, unsigned int height, unsigned int width) {
    printf("0 > Exporting rendering...\n");
    FILE *fp = fopen("rendering.txt","w");
    int i, j;
    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++)
            fprintf(fp, "%hhu ", rendering[i+j*height]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("0 > Export complete\n");
}

unsigned char computePixel(double complex d, unsigned int b, unsigned int N) {  //must specify complex
    unsigned char value = 1;
    long double complex z = 0 + 0*I;

    while((cabs(z) < b) && (value < N)) {
        z = z*z+d;
        ++value;
    }

    return value;
}

void computeRendering(unsigned int P, unsigned int N, unsigned int b, unsigned int height,
        unsigned int width, unsigned char * rendering, unsigned int p) {
    double dx = 2.0*b/(width-1);
    double dy = 2.0*b/(height-1);
    unsigned int wp = width/P;
    unsigned int hp = height;
    unsigned int xoff = p*width/P;
    unsigned int yoff = 0;
    unsigned int x, y, i = 0;
    double complex d;  // inbuilt function for dreal, dimag, 

    printf("%d > Processing rendering (%d:%d, %d:%d)...\n", p, xoff, xoff+wp, yoff, yoff+hp);

    for(x = 0; x < wp-1; ++x) {
        real(d) = (x+xoff)*dx-b;

        for(y = 0; y < hp-1; ++y) {
            imag(d) = (y+yoff)*dy-b;
            //d = dreal+I*dimag; inbuilt function for dreal and dimag

            rendering[y+x*hp] = computePixel(d, b, N);
        }
    }

    printf("%d > Rendering complete\n", p);
}

int main(int argc, char **argv) {
    unsigned int N = 196; //TODO: Final depth = 256
    unsigned int b = 2;
    unsigned int P;
    unsigned int height = 96; //TODO: final rendering in 2048*2048
    unsigned int width = 96;
    unsigned int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(width%P != 0) {
        printf("The number of processes (%d) must evenly divide the image width (%d)\n", P, width);
        return;
    }

    unsigned int size = height*width;
    unsigned int bufferSize = height*width/P;
    unsigned int wp = width/P;
    unsigned int hp = height;
    unsigned int offset;
    int dreal, dimag, d, q, i, j, x, y;

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
            printf("%d > Received rendering from %d\n", rank, q);
            offset = q*bufferSize;

            for(x = 0; x < wp-1; ++x) {
                for(y = 0; y < hp-1; ++y) {
                    rendering[y+x*hp+offset] = buffer[y+x*hp];
                }
            }
        }
        exportRendering(rendering, height, width);
    } else {
        MPI_Send(rendering, size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
