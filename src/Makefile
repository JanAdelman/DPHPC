# Since the blah target is first, it is the default target and will be run when we run "make"
compile:
	mpicxx -std=c++2a main_mpi.cpp -o main

run:
	mpiexec --oversubscribe -np 4 ./main --mca orte_base_help_aggregate 0

compile-run: 
	mpicxx -std=c++2a main_mpi.cpp -o main
	mpiexec --oversubscribe -np 4 ./main --mca orte_base_help_aggregate 0

compile-run-old: 
	mpicxx -std=c++11 main_mpi_cpy.cpp -o main
	mpiexec --oversubscribe -np 2 ./main --mca orte_base_help_aggregate 0

debug:
	mpicxx -g -std=c++2a main_mpi.cpp -o main
	mpiexec -n 4 gdb ./main

clean:
	rm -f main
