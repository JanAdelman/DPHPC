# Since the blah target is first, it is the default target and will be run when we run "make"
compile:
	mpicxx -std=c++2a -stdlib=libc++ main_mpi.cpp -o main

run:
	mpiexec -np 4 ./main --mca orte_base_help_aggregate 0

compile-run: 
	mpicxx -std=c++2a -stdlib=libc++ main_mpi.cpp -o main
	mpiexec -np 4 ./main --mca orte_base_help_aggregate 0

clean:
	rm -f main