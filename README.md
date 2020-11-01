# DPHPC

In this project we will implement a parallelised suffix tree building algorithm for fast string search (FSS). Our practical intentions are fast genome indexing. 

Inspiration is drawn from the following papers:

Aziz M.M.A., Thulasiraman P., Mohammed N. (2020) Parallel Generalized Suffix Tree Construction for Genomic Data. In: Martín-Vide C., Vega-Rodríguez M., Wheeler T. (eds) Algorithms for Computational Biology. AlCoB 2020. Lecture Notes in Computer Science, vol 12099. Springer, Cham. https://doi.org/10.1007/978-3-030-42266-0_1

Tian, Y., Tata, S., Hankins, R.A. et al. Practical methods for constructing suffix trees. The VLDB Journal 14, 281–299 (2005). https://doi.org/10.1007/s00778-005-0154-8

Our project consists of the following steps:

- Implement suffix tree building algorithm
- Implement merging of suffix tree 
- We will use MPI to coordinate across multiple working nodes

Remark:

In a shared memory system with 8 cores using vertical splitting (splitting the sequences vertically into 8 subsequences if 8 eight cores are used) the group of Al Aziz achieved a speedup of 4.45 compared to the sequential implementation. This is an indication that parallelising fast string search by employing suffix trees works. 

## Compilation

The build environment with MPICH and google-benchmark can be built using the contained Dockerfile. 

```
docker build -t hpc .
sudo docker run -it --entrypoint /bin/bash <ID>
```

The project is then built using MPI with the appropriate flags. 

```
mpicxx naive_suffixarray.cpp -o main -lbenchmark
mpiexec -np 1 ./main
```