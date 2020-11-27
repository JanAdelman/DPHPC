#include <iostream>
#include <string>
#include <mpi.h>
#include <math.h>
#include "helper.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

#define MASTER 0

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    //SETUP TUPLE STRUCT
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K};
    const MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //MAIN LOGIC

    char input[] = "agaagccagtactgcgacaaaggtaggacatggcgttgcaccaaatcagtaccggctccacaataattacaccatagggcaccgctatccgcgtgcgtca$";
    //char input[] = "ABCDEFGHAHAHABCDEFGHAHAHABCDEFGHAHAH$";
    int step_size = ceil((float)strlen(input) / (float)(world_size));//CHANGED TO EXLUCDE MASTER

    /*
    struct tuple_t buffer1;
    buffer1.idx = 0;
    strncpy(buffer1.seq, "ZEE", K);
    //buffer1.seq[K] = '\0';

    struct tuple_t buffer2;
    buffer2.idx = 1;
    strncpy(buffer2.seq, "EYE", K);
    //buffer2.seq[K] = '\0';

    tuple_t data[2] = {buffer1, buffer2};

    tuple_t_sort(data, 2);

    tuple_t_print(data, 2);
    */

    if (world_rank == MASTER)
    {

        //Send out overlapping substrings to workers

        //std::cout<<"world_rank"<<world_rank<<std::endl;

        int reciever = 1;
        //std::cout<<"thing"<<strlen(input) - step_size<<std::endl;
        for (int i = 0; i < strlen(input); i+=step_size)
        {
            if (i + step_size >= strlen(input)) //maybe add -1 oder >=
            { // Send rest
                int size = strlen(input)-i-K+1;

                tuple_t kmers[size];
                get_kmers(&input[i], K, kmers, size);
                tuple_t_sort(kmers, size);
                tuple_t* global_result_kmers;
                global_result_kmers = typename_t_sort<tuple_t>(log2(world_size), world_rank, kmers,size , MPI_COMM_WORLD);

                tuple_t_print(global_result_kmers, strlen(input)-K+1);

            }
            else
            {
                MPI_Send(&input[i], step_size + K - 1, MPI_CHAR, reciever, 0, MPI_COMM_WORLD);
            }
            reciever++;
        }
    }
    else
    {
        char sub_input[step_size + K - 1];

        int recieved_size;
        MPI_Status status;
        MPI_Recv(&sub_input, step_size + K - 1, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &recieved_size);

        tuple_t kmers[recieved_size-K+1];
        get_kmers(sub_input, K, kmers, recieved_size-K+1);

        // sort the tuples lexicographically
        tuple_t_sort(kmers, recieved_size-K+1);

        //Global sort
        typename_t_sort<tuple_t>(log2(world_size), world_rank, kmers, recieved_size-K+1, MPI_COMM_WORLD);

    }
    // Finalize the MPI environment.
    MPI_Finalize();
}
