CC = icc 

SOURCES = fft.c
EXEC = fft

CFLAGS = -I$(FFTW_HOME)/include -O3
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3_mpi -lfftw3
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

$(EXEC): $(SOURCES)
	$(CC) $(MPI_COMPILE_FLAGS) $(CFLAGS) $(SOURCES) $(MPI_LINK_FLAGS) $(LDFLAGS) -o $(EXEC)
clean:
	rm $(EXEC)
