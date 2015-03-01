/* Reaction-diffusion equation in 1D
 * Solution by Jacobi iteration
 * simple MPI implementation
 *
 * C Michael Hanke 2006-12-12
 */


// #define MIN(a,b) ((a) < (b) ? (a) : (b)) Not needed 

/* Use MPI */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
//Include math.h ??? 
#include "math.h"

/* define problem to be solved */
#define N 100   /* number of inner grid points */
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
//I will force this to be zero ie P will evenly divide N
    //I = (N+P-p-1)/P;   // (number of local elements) HOW CAN THIS BE, IT IS NOT EVEN AN INTEGER, HANKE WHAT HAVE THOUST DONE?
   
   // double n = p*L+MIN(p,R)+i; //(global index for given (p,i) COMPLICATION WHY?
    
/* Initialize MPI */
    unsigned int P;
    unsigned int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    int  L = N/P;
    int I = L; 
    double  h = 1/N;
    double R = N%P;  
    
    if (N < P) {
	fprintf(stdout, "Too few discretization points...\n");
	exit(1) ; 
    }

/* Compute local indices for data distribution, I DO NOT UNDERSTAND THIS (YET)*/
//double localI = p*L unecessary?
/* arrays */	
    double * unew;
    
    unew = (double *) malloc(I*sizeof(double));
/* Note: The following allocation includes additionally:
   - boundary conditins are set to zero
   - the initial guess is set to zero */
    double * u;  
    u = (double *) calloc(I+2, sizeof(double));

double ff[I]; // why does he suggest that this be external? I DO not understand...
double rr[I];
int 	k;
	for(k=1;k<I+1;k++){
		ff[k] = exp(h*(k+p*I)-0.5);
		rr[k] = sin(h*(k+p*I));
	}

/* Jacobi iteration */
    int step;
for (step = 0; step < SMX; step++) {
/* RB communication of overlap */
	if(p==0){ // First is red
		MPI_Send( (void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );
		printf("Process %d sent to process %d \n", p,  p +1);
		MPI_Recv( (&(u[I-1])),  1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
		printf("Process %d received from process %d \n", p,  p +1);

		u[0] = 0;
	}

	else if(p < P-1 && p %2!=0){ //Black blocks
		
		MPI_Recv(( &(u[0])),  1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
		printf("Process %d received from process %d \n", p,  p -1 );
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD);	
		printf("Process %d sent to process %d \n", p,  p -1);
		MPI_Recv(( &(u[I-1]) ), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE); // undefined for last process? whatevah
		printf("Process %d received from process %d \n", p,  p +1);
		MPI_Send((void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD);
		printf("Process %d sent to process %d \n", p,  p +1);
	
	}
	else if(p % 2 == 0 && p != 0){ // RED
		MPI_Send( (void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );
	    printf("Process %d sent to process %d \n", p,  p +1 );
		MPI_Recv( (&(u[I-1])),  1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );	
	printf("Process %d received from process %d \n", p,  p +1);
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );	
		printf("Process %d sent to process %d \n", p,  p -1);
	
		MPI_Recv((&(u[0])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
	}
	else{
		MPI_Recv(( &(u[0])),  1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
			printf("Process %d received from process %d \n", p,  p -1 );
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD);
		u[I-1]=0;
	
		
	}
	
/* local iteration step */
	int i;
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


unsigned int token;
if (p != 0) {
    MPI_Recv(&token, 1, MPI_INT, p - 1, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process %d received token %d from process %d\n",
           p, token, p - 1);
  FILE *fp = fopen("data.txt","r+");
  int j=0;
  for(j=0;j<I;j++){
  	fprintf(fp, "%hhu ", u[j+1]);
  	}
  fprintf(fp, "\n");
}
else {
    // Set the token's value if you are process 0
    token = 1;
    printf("Process %d set token to %d \n",
           p, token);
    printf(" Exporting data...\n");
    FILE *fp = fopen("data.txt","w");
    int j=0;	
    for(j=0;j<I;j++){
	fprintf(fp, "%hhu ", u[j+1]);
	}
  fprintf(fp, "\n");
}
MPI_Send(&token, 1, MPI_INT, (p + 1) % P,
         0, MPI_COMM_WORLD);

/* That's it */
    MPI_Finalize();
    exit(0);
}
