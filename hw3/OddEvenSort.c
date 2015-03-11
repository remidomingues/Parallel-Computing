
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

double cmpfunc (const void * a, const void * b)
{
   return ( *(double*)a - *(double*)b );
}



int main(int argc, char **argv)
{
unsigned int P;
unsigned int myrank;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &P);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

/* Find problem size N from command line */ //Remi could you please do this? shouldn't take you more than a second 
int N = 10;

//if (argc < 2) error_exit(“No size N given”);
//N = atoi(argv[1]);

/* local size. Modify if P does not divide N */

int R = N%P;
int I = (N-R)/P; //Number of local elements common to all
if( (R-myrank)>0 ){
  I = I+1;
  }


/* random number generator initialization */
srandom(myrank+1);
/* data generation */
double x[I];
int i;
for (i = 0; i < I; i++){
  x[i] = ((double) random())/(RAND_MAX+1);
}

// Local sort
qsort(x, I, sizeof(double), cmpfunc);

// Now we need to compare...

bool evenprocess = (myrank%2 == 0);
bool evenphase = true;
int step;

for(step=0; step<P;step++){  //As many steps as processors, though I'm not to sure about this...

  bool done = false;

  while ( done == false) {

    if( (evenphase  && evenprocess  )  xor ( !evenphase  && !evenprocess ) && myrank< P-1 ){
      
      MPI_Recv(&X, 1, MPI_DOUBLE, p + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Receives lowest from above neighbour
      //Now we need to check if this element has a place on this processor (is it smaller than the some element?)  
        bool b = false;
        int j = k; 
        double temp = X;
        
        while( b == false && j<I-1){ 
          if(X<x[j]){
            temp = x[j];
            x[j] = X;
            b = true;
          }
          j++;
          }
        
        //Elements need be exchanged
        
        MPI_Send( (void *) ( &(temp) ), 1, MPI_DOUBLE, p+1,1, MPI_COMM_WORLD);
    
    else if(myrank>0){
      MPI_Send((void *) (&(x[0])), 1, MPI_DOUBLE, p-1,1, MPI_COMM_WORLD); //Always send the smallest element
      MPI_Recv(&X, 1, MPI_DOUBLE, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if(x[0] == X ){
        done = true;  
      }
      else{   //Add new element to list and sort it in one step using bubble sort
        x[0] = X; //Added to list
        //Here follows bubble sort
        int j = 0;
        bool issorted = false;
      
        while( issorted == false ){
          issorted = true
          if(x[j+1]>x[j]){
            double temp = x[j];
            x[j] = x[j+1];
            x[j+1] = temp;
            issorted == false;
          }
          j++;
        }
      }
      
      }
    }
    
   evenphase != evenphase;
    
}
  
  




}
