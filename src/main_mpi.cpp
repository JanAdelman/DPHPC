#include <iostream>
#include <string>
#include <mpi.h>

#include "helper.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

#define MASTER 0
#define K 5

struct tuple_t
{
    int idx;
    char seq[K];
};

int main(int argc, char** argv){

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);



    MPI_Datatype tuple_struct;
    int lengths[2] = { 1, K };
    const MPI_Aint displacements[2] = { 0, sizeof(int)};
    MPI_Datatype types[2] = { MPI_INT, MPI_CHAR };
    MPI_Type_create_struct(2, lengths, displacements, types, &tuple_struct);
    MPI_Type_commit(&tuple_struct);


    //MAIN LOGIC

    if (world_rank == MASTER){
        
        struct tuple_t buffer1;
        buffer1.idx = 0;
        strncpy(buffer1.seq, "ADA", K);
        buffer1.seq[K] = '\0';

        struct tuple_t buffer2;
        buffer2.idx = 1;
        strncpy(buffer2.seq, "TZA", K);
        buffer2.seq[K] = '\0';

        tuple_t data[2] = {buffer1, buffer2}; 

        MPI_Send(&data, 2, tuple_struct, 1, 0, MPI_COMM_WORLD);
    
    }
    if (world_rank == 1){
        
        tuple_t r_data[2]; 

        MPI_Recv(&r_data, 2, tuple_struct, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

        std::cout << r_data[0].idx << r_data[0].seq << std::endl;
    
    }

    //scattering str from MASTER node out to the other node
	//MPI_Scatter(b, n_per_proc, MPI_INT, bp, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);

    //compute on each substring
    //tuple_vector sub_vec =  get_kmers(substr, K);

    //MASTER node gathering array from the workers
    //MPI_Gather(cp, n_per_proc, MPI_INT, c, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);


    // Finalize the MPI environment.
    MPI_Finalize();

}