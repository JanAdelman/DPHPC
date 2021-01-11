#!/bin/bash

mpicxx -std=c++11 -O3 main_mpi.cpp -o main

echo "" > performance.csv
echo "String_lenght, #Cores, Repetition, Overall, Min, Max, Average" > performance.csv 

declare -a cores=(32 64 128)

for file in ./data/*
do
for i in ${cores[@]}
do
for numbers in {1..10}
do
arr_settings+=($file $i $numbers)
echo -e "${arr_settings[0]},${arr_settings[1]}, ${arr_settings[2]},  \c" >> performance.csv 
mpirun -n ${i} ./main ${file} >> performance.csv 
declare -a arr_settings=()
done
done
done
