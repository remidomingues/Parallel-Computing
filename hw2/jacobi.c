
/* Use MPI */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

/* define problem to be solved */
#define N 1000   /* number of inner grid points */
#define SMX 1000000 /* number of iterations */


int main(int argc, char *argv[])
{

/* Initialize MPI */
    unsigned int P;
    unsigned int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    int  L = N/P;
    int I = L; 
    double   h = 1.0/(N+1); // We want N inner points
    double R = N%P;  
    
    if (N < P) {
	fprintf(stdout, "Too few discretization points...\n");
	exit(1) ; 
    }
/* arrays */	
    double * unew;
    
    unew = (double *) malloc(I*sizeof(double));

    double * u;  
    u = (double *) calloc(I+2, sizeof(double));

double ff[I];
double rr[I];
int 	k;
	for(k=0;k<I;k++){
		ff[k] = -1/pow(h*(k+1+p*I),1.5);
		rr[k] = pow(h*(k+1+p*I),3);
	}
/* Jacobi iteration */
    int step;
for (step = 0; step < SMX; step++) {
/* RB communication of overlap */
	if(p==0){ // First is red and contains left boundary 
		MPI_Send( (void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );
		MPI_Recv( (&(u[I-1])),  1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
		u[0] = 0;
	}

	else if(p < P-1 && p %2!=0){ //Black blocks
		
		MPI_Recv(( &(u[0])),  1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD);	
		MPI_Recv(( &(u[I-1]) ), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
		MPI_Send((void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD);
	
	}
	else if(p % 2 == 0 && p != 0){ // Red blocks
		MPI_Send( (void *) (&(u[I-2])), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD );
	   	MPI_Recv( (&(u[I-1])),  1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );	
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD );	
		MPI_Recv((&(u[0])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
	}
	else{  //Last block contains right boundary condition
		MPI_Recv(( &(u[0])),  1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
		MPI_Send((void *) (&(u[1])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD);
		u[I-1]=0;
	
	}
	
/* local iteration step */
	int i;
	for(i=0; i<I;i++ ){
	    unew[i] = (u[i]+u[i+2]-h*h*ff[i])/(2.0-h*h*rr[i]);
	}
	for(i=0;i<I;i++){  
	    u[i+1] = unew[i];
	    }
}

unsigned int token;
if (p != 0) {
    MPI_Recv(&token, 1, MPI_INT, p - 1, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process %d received token %d from process %d\n",p, token, p - 1);
  FILE *fp = fopen("data.txt","a");
  int j=0;
  for(j=0;j<I;j++){
  	fprintf(fp, "%f ", u[j+1]);
  	}
  fprintf(fp, "\n");
  fclose(fp);
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
	fprintf(fp, "%f ", u[j+1]);
	}
  fprintf(fp, "\n");
  fclose(fp);
}
MPI_Send(&token, 1, MPI_INT, (p + 1) % P,
         0, MPI_COMM_WORLD);

/* That's it */
    MPI_Finalize();
    exit(0);
}
