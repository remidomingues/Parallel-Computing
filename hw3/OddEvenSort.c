int main(int argc, char **argv)
{
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &P);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
/* Find problem size N from command line */
if (argc < 2) error_exit(“No size N given”);
N = atoi(argv[1]);
/* local size. Modify if P does not divide N */
I = N/P;
/* random number generator initialization */
srandom(myrank+1);
/* data generation */
for (i = 0; i < I; i++)
  x[i] = ((double) random())/(RAND_MAX+1);

}
