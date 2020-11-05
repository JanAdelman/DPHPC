# DPHPC

In this project we will implement a parallelised suffix array building algorithm for fast string search (FSS). Our practical intentions are fast genome indexing. 

Inspiration is drawn from the following papers:

Flick, Patrick, and Srinivas Aluru. “Parallel Distributed Memory Construction of Suffix and Longest Common Prefix Arrays.” Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis on - SC ’15, ACM Press, 2015, pp. 1–10. DOI.org (Crossref), https://doi.org/10.1145/2807591.2807609.

Aziz M.M.A., Thulasiraman P., Mohammed N. (2020) Parallel Generalized Suffix Tree Construction for Genomic Data. In: Martín-Vide C., Vega-Rodríguez M., Wheeler T. (eds) Algorithms for Computational Biology. AlCoB 2020. Lecture Notes in Computer Science, vol 12099. Springer, Cham. https://doi.org/10.1007/978-3-030-42266-0_1

Tian, Y., Tata, S., Hankins, R.A. et al. Practical methods for constructing suffix trees. The VLDB Journal 14, 281–299 (2005). https://doi.org/10.1007/s00778-005-0154-8

Our project consists of the following steps:

- Implement suffix array building algorithm
- We will use MPI to coordinate across multiple working nodes

## Compilation

The build environment with MPICH and google-benchmark can be built using the contained Dockerfile. 

```
docker build -t hpc .
sudo docker run -it --entrypoint /bin/bash <ID>
```

The project is then built using MPI with the appropriate flags. 

```
mpicxx -std=c++14 -stdlib=libc++ prefix_doubling.cpp -o main -lbenchmark
mpiexec -np 1 ./main
```
