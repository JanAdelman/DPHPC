# Since the blah target is first, it is the default target and will be run when we run "make"
compile:
	mpicxx -std=c++2a main_mpi.cpp -o main

run:
	mpiexec --oversubscribe -np 4 ./main --mca orte_base_help_aggregate 0

compile-run: 
	mpicxx -std=c++11 -O3 main_mpi.cpp -o main
	mpiexec -np 4 ./main

debug:
	mpicxx -g -std=c++2a main_mpi.cpp -o main
	mpiexec -n 4 gdb ./main

clean:
	rm -f main
