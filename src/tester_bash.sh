#!/bin/bash

mpicxx -std=c++11 -O3 main_mpi.cpp -o main
echo "" > performance.csv

#cores=(32,64,128,256,512,1024)
declare -a cores=(2 4)
for file in ./data/*
do
for i in ${cores[@]}
do
for numbers in {1..2}
do
arr_settings+=($file $i $numbers)
echo "${arr_settings[0]},${arr_settings[1]}, ${arr_settings[2]},  \c" >> performance.csv
mpiexec -np $i ./main $file >> performance.csv 
declare -a arr_settings=()
done
done
done
