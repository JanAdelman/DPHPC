#include <mpi.h>
#include <stdio.h>

#include "helper.cpp"

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //int sendbuf[] = {0,0,2,2,2,5,5,5,8,8};
    int sendbuf[] = {0,1,2,3,4,5,6,7,8,9};
    int h = 6; 
    int* recv_data; 
    
    int local_length;

    if (world_rank == 0){
         local_length = 3;
         recv_data = (int *) malloc(sizeof(int) * local_length);

         for (int i = 0; i<local_length; i++)
            recv_data[i] = sendbuf[0+i];
        ///print_int_array(recv_data, 4);
    }
    if (world_rank == 1){
        local_length = 3;
        recv_data = (int *) malloc(sizeof(int) * local_length);
         for (int i = 0; i<local_length; i++)
            recv_data[i] = sendbuf[3+i];
        }
    if (world_rank == 2){
        local_length = 2;
        recv_data = (int *) malloc(sizeof(int) * local_length);
                 for (int i = 0; i<local_length; i++)
            recv_data[i] = sendbuf[6+i];
     }
        
    if (world_rank == 3){
        local_length = 2;
        recv_data = (int *) malloc(sizeof(int) * local_length);
                 for (int i = 0; i<local_length; i++)
            recv_data[i] = sendbuf[8+i];
    }
    
    

    int new_idx[4];
    int offsets[4] = {2,5,7,9};
    int s_offsets[4] = {0,3,6,8};

    naive_shift(recv_data, h, MPI_COMM_WORLD, world_rank,
            world_size,s_offsets , local_length, offsets);

    // Finalize the MPI environment.
    MPI_Finalize();
}
