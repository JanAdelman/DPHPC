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

    //char input[] = "agaagccagtactgcgacaaaggtaggacatggcgttgcaccaaatcagtaccggctccacaataattacaccatagggcaccgctatccgcgtgcgtca$";
    //char input[] = "ABCDEFGHAHAHABCDEFGHAHAHABCDEFGHAHAH$";
    char input[] = "AAAASBBBL"; //overlapping string of k=2 and p=4
    //AAA AS BB BL
    //AAAA ASB BBB BL
    //char input[]="AAAAIBBBBCDEFGHI";
    int string_length = strlen(input);
    //int step_size = ceil((float)string_length / (float)(world_size)); //CHANGED TO EXLUCDE MASTER
    int nmin = string_length / world_size;   //min step size to send
    int nextra = string_length % world_size; //remainder if not divisible by nprocessors
    //std::cout<<nmin<<std::endl;
    //std::cout<<nextra<<std::endl;
    int k = 0;
    int sendnums[world_size]; //size to send out to process i
    int displace[world_size];  //displacement from start
    //maybe do this world_size-1 to not include 0
    for (int i = 1; i < world_size; i++)
    {
        if (i - 1 < nextra) sendnums[i] = nmin + 1;
        else sendnums[i] = nmin;
        displace[i] = k;
        k = k + sendnums[i];
    }
    sendnums[0] = nmin;
    displace[0] = k;
    tuple_t *global_result_kmers;

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

        
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Send(&input[displace[rank]], sendnums[rank] + K - 1, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        }
        int size = sendnums[0] - K + 1;
        tuple_t kmers[size];
        get_kmers(&input[displace[0]], K, kmers, size);
        tuple_t_sort(kmers, size);
        //tuple_t_print(kmers, size);
        global_result_kmers = typename_t_sort<tuple_t>(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
        tuple_t_print(global_result_kmers, string_length - K + 1);

        /*
        int reciever = 1;
        //std::cout<<"thing"<<strlen(input) - step_size<<std::endl;
        for (int i = 0; i < string_length; i += step_size)
        {
            std::cout<<i + step_size<<std::endl;
            if (i + step_size >string_length-1) //maybe add -1 oder >=
            {                                   // Send rest
            std::cout<<"if statement"<<std::endl;
                int size = string_length - i - K + 1;
                tuple_t kmers[size];
                get_kmers(&input[i], K, kmers, size);
                tuple_t_sort(kmers, size);
                std::cout<<"kmers sorted"<<std::endl;
                global_result_kmers = typename_t_sort<tuple_t>(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
                std::cout<<"sort this shit"<<std::endl;
                tuple_t_print(global_result_kmers, string_length-K+1);
            }
            else
            {
                MPI_Send(&input[i], step_size + K - 1, MPI_CHAR, reciever, 0, MPI_COMM_WORLD);
            }
            reciever++;
        }
        */
    }
    else
    {
        int size = sendnums[world_rank]+K-1;
        //std::cout<<world_rank<<std::endl;
        //std::cout<<"size "<<size<<std::endl;
        char sub_input[size];

        //int recieved_size;
        //MPI_Status status;
        MPI_Recv(&sub_input, size,MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, NULL);
        //MPI_Get_count(&status, MPI_CHAR, &recieved_size);

        tuple_t kmers[size - K + 1];
        get_kmers(sub_input, K, kmers, size - K + 1);
        //std::cout<<world_rank<<std::endl;
        //tuple_t_print(kmers, size - K + 1);

        // sort the tuples lexicographically
        tuple_t_sort(kmers, size - K + 1);
        //std::cout<<world_rank<<std::endl;
        //tuple_t_print(kmers, size - K + 1);

        //Global sort
        typename_t_sort<tuple_t>(log2(world_size), world_rank, kmers, size - K + 1, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int kmin = (string_length - K + 1) / world_size;   //min step size to send
    int kextra = (string_length - K + 1) % world_size; //remainder if not divisible by nprocessors
    int now = 0;
    int sendcounts[world_size]; //size to send out to process i
    int displs[world_size];     //displacement from start
    //maybe do this world_size-1 to not include 0
    for (int i = 0; i < world_size; i++)
    {
        if (i < kextra)
            sendcounts[i] = kmin + 1;
        else
            sendcounts[i] = kmin;
        displs[i] = now;
        now = now + sendcounts[i];
    }
    tuple_t recvbuf[sendcounts[world_rank]];
    MPI_Scatterv(global_result_kmers, sendcounts, displs, MPI_TUPLE_STRUCT, recvbuf, sendcounts[world_rank], MPI_TUPLE_STRUCT, 0, MPI_COMM_WORLD);
    int local_array[sendcounts[world_rank]];
    int local_SA[sendcounts[world_rank]];
    rebucketing(local_array, recvbuf, sendcounts[world_rank], displs[world_rank], local_SA);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank < world_size - 1)
    {
        tuple_t final[1];
        final[0].idx = local_array[sendcounts[world_rank] - 1];
        memcpy(final[0].seq, recvbuf[sendcounts[world_rank] - 1].seq, sizeof(recvbuf[0].seq));
        //final[0].seq=recvbuf[sendcounts[world_rank]-1].seq;
        MPI_Send(&final, 1, MPI_TUPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        tuple_t curr[1];
        MPI_Recv(&curr, 1, MPI_TUPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        //tuple_t_print(curr, K);
        int i = 0;
        //std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        //std::cout<<curr[0].idx<<std::endl;
        while (char_array_comp(curr[0].seq, recvbuf[i].seq, K))
        {
            //std::cout<<"while_loop"<<std::endl;
            local_array[i] = curr[0].idx;
            i++;
        }
    }
    //free(recvbuf);
    MPI_Barrier(MPI_COMM_WORLD);
    /*
    for (int i = 0; i < sendcounts[world_rank]; i++)
    {
        std::cout << local_SA[i] << std::endl;
    }
    */
    /*
    int global_array[string_length-K+1];
    MPI_Gatherv(local_array,sendcounts[world_rank], MPI_INT,global_array, sendcounts,displs,MPI_INT,0,MPI_COMM_WORLD);
    
    if(world_rank==0){
        for(int i=0;i<string_length-K+1;i++){
        std::cout<<global_array[i]<<",";
    }
    }
    */

    // Finalize the MPI environment.
    MPI_Finalize();
}
