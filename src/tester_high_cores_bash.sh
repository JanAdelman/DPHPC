#!/bin/bash

mpicxx -std=c++11 -O3 main_mpi.cpp -o main
echo "" > big_inputs.csv
echo "String_length, #Cores, Repetition, Overall, Min, Max, Average " > big_inputs.csv 

declare -a cores=(128 256 512 1024 2048)
for file in ./data_big/*
do
for i in ${cores[@]}
do
for numbers in {1..10}
do
arr_settings+=($file $i $numbers)
echo -e "${arr_settings[0]},${arr_settings[1]}, ${arr_settings[2]}, \c" >> big_inputs.csv
mpiexec -n $i ./main $file >> big_inputs.csv 
declare -a arr_settings=()
done
done
done
