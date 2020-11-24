#include <iostream>
#include <string>
#include <mpi.h>

#include "helper.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

#define MASTER 0

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Datatype tuple_struct;
    int lengths[2] = {1, K};
    const MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &tuple_struct);
    MPI_Type_commit(&tuple_struct);

    //MAIN LOGIC

    char input[] = "ABCDEFGHIJKLMNOPQRST$";
    int step_size = strlen(input) / (world_size - 1);//CHANGED TO EXLUCDE MASTER
    //std::cout<<"step: "<<step_size<<std::endl;

    if (world_rank == MASTER)
    {

        
        /*
        struct tuple_t buffer1;
        buffer1.idx = 0;
        strncpy(buffer1.seq, "ADA", K);
        buffer1.seq[K] = '\0';

        struct tuple_t buffer2;
        buffer2.idx = 1;
        strncpy(buffer2.seq, "TZA", K);
        buffer2.seq[K] = '\0';

        tuple_t data[2] = {buffer1, buffer2}; 
        */

        //Send out overlapping substrings to workers

        //std::cout<<"world_rank"<<world_rank<<std::endl;

        int reciever = 1;
        //std::cout<<"thing"<<strlen(input) - step_size<<std::endl;
        for (int i = 0; i <= strlen(input) - step_size; i+=step_size)
        {
            //std::cout<<"receiver "<<reciever<<std::endl;
            if (i + step_size + K > strlen(input))
            { // Send rest
            /*
            std::cout<<i+step_size+K<<std::endl;
            std::cout<<strlen(input)<<std::endl;
            std::cout<<"if "<<std::endl;
            std::cout<<i<<std::endl;
            std::cout<<strlen(input)-i<<std::endl;
            */
                MPI_Send(&input[i], strlen(input) - i, MPI_CHAR, reciever, 0, MPI_COMM_WORLD);
            }
            else
            {
                //std::cout<<"else "<<std::endl;
                MPI_Send(&input[i], step_size + K - 1, MPI_CHAR, reciever, 0, MPI_COMM_WORLD);
            }
            reciever++;
        }
    }
    else
    {
        //std::cout<<"world_rank"<<world_rank<<std::endl;
        char sub_input[step_size + K - 1];

        int recieved_size;  
        MPI_Status status;
        MPI_Recv(&sub_input, step_size + K - 1, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &recieved_size);

        //std::cout << "WHOLE SUB INPUT: "; print_char_array(sub_input);

        tuple_t kmers[recieved_size-K+1];
        get_kmers(sub_input, K, kmers, recieved_size-K+1);

        //MPI_Barrier(MPI_COMM_WORLD); 
        print_kmers(kmers, recieved_size-K+1);

    }
    // Finalize the MPI environment.
    MPI_Finalize();
}