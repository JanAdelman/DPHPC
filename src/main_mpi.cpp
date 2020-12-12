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

    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    const MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
    MPI_Datatype types_ISA[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths_ISA, displacements_ISA, types_ISA, &MPI_TUPLE_ISA);
    MPI_Type_commit(&MPI_TUPLE_ISA);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //MAIN LOGIC

    //char input[] = "agaagccagtactgcgacaaaggtaggacatggcgttgcaccaaatcagtaccggctccacaataattacaccatagggcaccgctatccgcgtgcgtca$";
    //char input[] = "ABCDEFGHAHAHABCDEFGHAHAHABCDEFGHAHAH$";
    //char input[] = "AAAASBBBL"; //overlapping string of k=2 and p=4
    //char input[]="NJSDJSNMCXNCJDKHFUIAAAJ";
    char input[]="MISSISSIPPI$";
    //AAA AAA AAH AHH HHH HHK HKK KKK KKK
    //0 1 2 3 4 5 6 7 8
    //0 1 2 3 
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
        get_kmers_adapt(&input[displace[0]], K, kmers, size,displace[0]);
        t_sort(kmers, size);
        //tuple_t_print(kmers, size);
        global_result_kmers = typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
        //t_print(global_result_kmers, string_length - K + 1);

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
        get_kmers_adapt(sub_input, K, kmers, size - K + 1,displace[world_rank]);
        //std::cout<<world_rank<<std::endl;
        //tuple_t_print(kmers, size - K + 1);

        // sort the tuples lexicographically
        t_sort(kmers, size - K + 1);
        //std::cout<<world_rank<<std::endl;
        //tuple_t_print(kmers, size - K + 1);

        //Global sort
        typename_t_sort(log2(world_size), world_rank, kmers, size - K + 1, MPI_COMM_WORLD);
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
    //std::cout<<"before scattering"<<std::endl;
    MPI_Scatterv(global_result_kmers, sendcounts, displs, MPI_TUPLE_STRUCT, recvbuf, sendcounts[world_rank], MPI_TUPLE_STRUCT, 0, MPI_COMM_WORLD);
    tuple_ISA SA_B[sendcounts[world_rank]];//array of tuples on this processor

    int counts_bucket[world_size];//counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    
    rebucketing(SA_B, recvbuf, sendcounts[world_rank], displs,counts_bucket,world_rank, world_size);
    //std::cout<<"rebucketing done"<<std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank < world_size - 1)
    {
        tuple_t final[1];
        final[0].idx = SA_B[sendcounts[world_rank] - 1].B;
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
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //########################################################
    //########################################################
    //########################################################
    //########################################################


    for(int h=K;h<=string_length-K+1;h*=2)
    {
    if(world_rank==0){
        std::cout<<"h"<<h<<std::endl;
    }


    
    tuple_ISA sendbufISA[sendcounts[world_rank]];
    int displace_ISA[world_size];
    int current_index_ISA=0;
    int counts_filled[world_size];//used to fill the sendbuf ISA
    for(int i=0;i<world_size;i++){
        displace_ISA[i]=current_index_ISA;
        current_index_ISA+=counts_bucket[i];
        counts_filled[i]=0;
    }

    //std::cout << "world rank after sendbuf made"<<world_rank << std::endl;

    for(int j=0;j<sendcounts[world_rank];j++){
        //int id= bucket_id(displs, SA_B[j], world_size);
        int id= bucket_id(displs, SA_B[j], world_size);
        sendbufISA[displace_ISA[id]+counts_filled[id]]=SA_B[j];
        counts_filled[id]++;    
    }

  


    //this array must be updated after every for loop
    tuple_ISA* recvbuf_ISA= probing_alltoallv(sendbufISA, displace_ISA, sendcounts[world_rank], world_size, counts_bucket, MPI_COMM_WORLD, world_rank, MPI_TUPLE_ISA);
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    
    int B[sendcounts[world_rank]];
    reorder_to_stringorder(B,recvbuf_ISA,sendcounts[world_rank],displs[world_rank]);

    

    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displs[i]+sendcounts[i]-1;
    }

    int B2[sendcounts[world_rank]];
    std::copy(B,B+sendcounts[world_rank],B2);

    //SHIFTING 
    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displs, sendcounts[world_rank], offsets);

    
    //TRIPLES 

    triple_t triple_arr[sendcounts[world_rank]];
    create_triple(B,B2,sendcounts[world_rank],displs[world_rank], triple_arr);
    //std::cout << "-----" << std::endl;
    
    //SORT LOCALLY
    t_sort(triple_arr, sendcounts[world_rank]);

    //SORT GLOBALLY
    triple_t* global_result_triplet;

    if (world_rank == MASTER){
        global_result_triplet = (triple_t*) malloc(sizeof(triple_t) * (string_length - K + 1));
        global_result_triplet = typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
    }
    else
    {
        typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
    }

    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    const MPI_Aint displacements_triple[3] = {0, sizeof(int), 2*sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);



    
    //SCATTER

    triple_t recvbuf_triplets[sendcounts[world_rank]];
    MPI_Scatterv(global_result_triplet, sendcounts, displs, MPI_TRIPLE_STRUCT, recvbuf_triplets, sendcounts[world_rank], MPI_TRIPLE_STRUCT, 0, MPI_COMM_WORLD);

    //REBUCKET

    

    //int counts_bucket[world_size];//counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    //tuple_ISA SA_B[sendcounts[world_rank]];
    rebucketing(SA_B,recvbuf_triplets, sendcounts[world_rank], displs,world_rank,world_size,counts_bucket);

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank < world_size - 1)
    {
        triple_t final[1];
        final[0].b= recvbuf_triplets[sendcounts[world_rank] - 1].b;
        final[0].b2= recvbuf_triplets[sendcounts[world_rank] - 1].b2;
        final[0].idx= SA_B[sendcounts[world_rank] - 1].B;
        //final[0].seq=recvbuf[sendcounts[world_rank]-1].seq;
        MPI_Send(&final, 1, MPI_TRIPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        triple_t curr[1];
        MPI_Recv(&curr, 1, MPI_TRIPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        //tuple_t_print(curr, K);
        int i = 0;
        //std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        //std::cout<<curr[0].idx<<std::endl;
        while (curr[0].b == recvbuf_triplets[i].b and curr[0].b2 == recvbuf_triplets[i].b2)
        {
            //std::cout<<"while_loop"<<std::endl;
            SA_B->B = curr[0].idx;
            i++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

   if(h==4){
      std::cout<<"world_rank"<<world_rank<<"SAB CREATED in for loop"<<std::endl;
    tuple_print(SA_B,sendcounts[world_rank]); 
   }
   
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, sendcounts[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout<<"World_rank"<<world_rank<<" h "<<h<<" singleton "<<singleton<<std::endl;

    bool singleton_global;
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);

    if(world_rank==0){  
    std::cout<<"End of one iteration, h:"<<h<<" "<<singleton_global<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(singleton_global && world_rank==0){
        std::cout << singleton_global << std::endl;
        std::cout << "final";
        break;
    }

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    // Finalize the MPI environment.
    MPI_Finalize();
    
    return 0;
}
