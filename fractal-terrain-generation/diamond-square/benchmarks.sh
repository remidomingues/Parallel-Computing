#!/bin/sh
cd ~/Downloads/tmp

for iter in `seq 8 13`; do
    printf "time_$iter = ["
    for n in `seq 1 8`; do
        time1=`nice -20 mpirun -np $n ./diamond_square $iter terrain.in`
        time2=`nice -20 mpirun -np $n ./diamond_square $iter terrain.in`
        time3=`nice -20 mpirun -np $n ./diamond_square $iter terrain.in`
        r=$((( $time1 + $time2 + $time3 ) / 3))
        printf "%.3f; " $r
    done
    printf "]\n"
done
