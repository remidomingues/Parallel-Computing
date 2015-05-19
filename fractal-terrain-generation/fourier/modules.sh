#!/bin/bash

ic='i-compilers/11.1';
mpi='openmpi/1.5.4-intel-11';
fftw='fftw/3.3.3-intel-11-dp-openmpi-1.5';

module add $ic $mpi $fftw;
LD_LIBRARY_PATH=${FFTW_HOME}/lib:${LD_LIBRARY_PATH}
