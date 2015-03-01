/* Reaction-diffusion equation in 1D
 * Solution by Jacobi iteration
 * simple MPI implementation
 *
 * C Michael Hanke 2006-12-12
 */


// #define MIN(a,b) ((a) < (b) ? (a) : (b)) Not needed 

/* Use MPI */
#include "mpi.h"

//Include math.h ??? 
#include "math.h"

/* define problem to be solved */
#define N 1000   /* number of inner grid points */
#define SMX 2 // 1000000 /* number of iterations */

/* implement coefficient functions */
extern double r(const double x);
extern double f(const double x);

/* We assume linear data distribution. The formulae according to the lecture
   are:
      L = N/P;
      R = N%P;
      I = (N+P-p-1)/P;    (number of local elements)
      n = p*L+MIN(p,R)+i; (global index for given (p,i)
   Attention: We use a small trick for introducing the boundary conditions:
      - The first ghost point on p = 0 holds u(0);
      - the last ghost point on p = P-1 holds u(1).
   Hence, all local vectors hold I elements while u has I+2 elements.
*/

int main(int argc, char *argv[])
{
/* local variable */
   double  L = N/P;
   double  h = 1/N;
   double R = N%P; //I will force this to be zero ie P will evenly divide N
    //I = (N+P-p-1)/P;   // (number of local elements) HOW CAN THIS BE, IT IS NOT EVEN AN INTEGER, HANKE WHAT HAVE THOUST DONE?
   double I = L; 
   // double n = p*L+MIN(p,R)+i; //(global index for given (p,i) COMPLICATION WHY?
    
/* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    
 
    if (N < P) {
	fprintf(stdout, "Too few discretization points...\n");
	exit(1);
    }

/* Compute local indices for data distribution, I DO NOT UNDERSTAND THIS (YET)*/
//double localI = p*L unecessary?
/* arrays */
    unew = (double *) malloc(I*sizeof(double));
/* Note: The following allocation includes additionally:
   - boundary conditins are set to zero
   - the initial guess is set to zero */
    u = (double *) calloc(I+2, sizeof(double));


/* Jacobi iteration */
    for (step = 0; step < SMX; step++) {
/* RB communication of overlap */
	if(p % 2 == 0){ // red?  From slides, TO DO
		MPI_send(u[I-2], L, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );
		//send(u[Ip_2],p+1);
		MPI_recv(u[I-1], L, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );	
		//receive(u[Ip-1],p+1);
		MPI_send(u[1], L, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );	// undefined for first process?
		//send(u[1],p-1);
		MPI_recv(u[0], L, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );
		//receive(u[0],p-1);
	else{
		MPI_recv(u[0], L, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );
		MPI_send(u[1], L, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );	
		MPI_recv(u[I-1], L, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD ); // undefined for last process? whatevah
		MPI_send(u[I-2], L, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );	
		/* receive(u[0],p-1);
		send(u[1],p-1);
		receive(u[Ip-1],p+1);
		send(u[Ip-2],p+1);*/
	}
	}

/* local iteration step */
	for(i=0; i<I;i++ ){
	    unew[i] = (u[i]+u[i+2]-h*h*ff[i])/(2.0-h*h*rr[i]);
	}
	for(i=0;i<I;i++){  //This is ugly coding, but i don't know any better... :(
	    u[i+1] = unew[i];
	    }
}

/* output for graphical representation */
/* Instead of using gather (which may lead to excessive memory requirements
   on the master process) each process will write its own data portion. This
   introduces a sequentialization: the hard disk can only write (efficiently)
   sequentially. Therefore, we use the following strategy:
   1. The master process writes its portion. (file creation)
   2. The master sends a signal to process 1 to start writing.
   3. Process p waites for the signal from process p-1 to arrive.
   4. Process p writes its portion to disk. (append to file)
   5. process p sends the signal to process p+1 (if it exists).
*/
int MPI_wait(MPI_Request *request, MPI_Status *status)
 
int	right = (p + 1) % numprocs;
int	left = p - 1;
    if (left < 0){
        left = numprocs - 1;
    }

if(p == 0) {
	
	printf("0 > Exporting data...\n");
	FILE *fp = fopen("data.txt","w");
	int j
	for(j=0;j<I;j++){
	fprintf(fp, "%hhu ", u[j+1]);
	}
	fprintf(fp, "\n");
	MPI_send(1, 1, MPI_INT, right, 123, MPI_COMM_WORLD);
}
else{
	int ok;
	MPI_Irecv(ok, 1, MPI_INT, left, 123, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	FILE *fp = fopen("data.txt","r+");
	int j
	for(j=0;j<I;j++){
		fprintf(fp, "%hhu ", u[j+1]);
	}
	fprintf(fp, "\n");
	MPI_send(1, 1, MPI_INT, right, 123, MPI_COMM_WORLD);
}



int i, j;
for (j = 0; j < width; j++) {
for (i = 0; i < height; i++)
fprintf(fp, "%hhu ", rendering[i+j*height]);
fprintf(fp, "\n");
}
fclose(fp);
printf("0 > Export complete\n");


/* That's it */
    MPI_Finalize();
    exit(0);
}
