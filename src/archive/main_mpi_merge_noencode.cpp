#include <iostream>
#include <string>
#include <mpi.h>
#include <math.h>
#include "helper_merge_noencode.cpp"
#include <fstream>
#include <vector>
#include "mem-usage.h"

#define MASTER 0

int main(int argc, char **argv)
{
    //setting up timing variables
    double start, end;
    MPI_Init(NULL, NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    // Starting timer
    start = MPI_Wtime();

    //SETUP TUPLE STRUCT
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K_SIZE};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
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
    //reading in file to get the string length
    std::ifstream is;
    is.open("./10power8.txt", std::ios::binary );
    is.seekg (0, std::ios::end);
    int string_length = is.tellg();
    is.close();

    int nmin = string_length / world_size;   //min step size to send
    int nextra = string_length % world_size; //remainder if not divisible by nprocessors

    int k = 0;
    int sendnums[world_size]; //size to send out to process i
    int displace[world_size];  //displacement from start
    for (int i = 0; i < world_size; i++)
    {
        if (i < nextra) sendnums[i] = nmin + 1;
        else sendnums[i] = nmin;
        displace[i] = k;
        k = k + sendnums[i];
    }


    std::fstream File("./10power8.txt", std::ios::in | std::ios::out );
    File.seekg(displace[world_rank], std::ios::beg);
    char* input;
    int size;
    if(world_rank==world_size-1){
        input=new char[sendnums[world_rank]];
        File.read(input, sendnums[world_rank]);
        size=sendnums[world_rank]-K_SIZE+1;
    }
    else{
        input=new char[sendnums[world_rank]+K_SIZE-1];
        File.read(input, sendnums[world_rank]+K_SIZE-1);
        size=sendnums[world_rank];
    }
    File.close();

    //Declare vector for global result
    std::vector<tuple_t> global_result_kmers; 

    if (world_rank == MASTER)
    {

        //get kmers
        std::vector<tuple_t> kmers(size);
        get_kmers_adapt(input, K_SIZE, kmers, size,displace[0]);
        delete[] input;
        //sorting tuples locally
        t_sort(kmers, size);

        //global sort
        global_result_kmers.resize(string_length - K_SIZE + 1);
        global_result_kmers = typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
 
        std::vector<tuple_t>().swap(kmers);

    }
    else
    {
       
        std::vector<tuple_t> kmers(size);
        get_kmers_adapt(input, K_SIZE, kmers, size,displace[world_rank]);

        delete[] input;

        // sort the tuples locally
        t_sort(kmers, size);

        //Global sort
        typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);

        std::vector<tuple_t>().swap(kmers);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int kmin = (string_length - K_SIZE + 1) / world_size;   //min step size to send
    int kextra = (string_length - K_SIZE + 1) % world_size; //remainder if not divisible by nprocessors
    int now = 0;
    int sendcounts[world_size]; //size to send out to process i
    int displs[world_size];     //displacement from start
    
    for (int i = 0; i < world_size; i++)
    {
        if (i < kextra)
            sendcounts[i] = kmin + 1;
        else
            sendcounts[i] = kmin;
        displs[i] = now;
        now = now + sendcounts[i];
    }


    std::vector<tuple_t> recvbuf(sendcounts[world_rank]);
    MPI_Scatterv(&global_result_kmers[0], sendcounts, displs, MPI_TUPLE_STRUCT,
     &recvbuf[0], sendcounts[world_rank], MPI_TUPLE_STRUCT, 0, MPI_COMM_WORLD);

    std::vector<tuple_t>().swap(global_result_kmers);

    std::vector<tuple_ISA> SA_B(sendcounts[world_rank]);

    int counts_bucket[world_size];//counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    //Rebucketing
    rebucketing(SA_B, recvbuf, sendcounts[world_rank], displs,counts_bucket,world_rank, world_size);
    MPI_Barrier(MPI_COMM_WORLD);


    if (world_rank < world_size - 1)
    {
        //sending last element to neighbor
        tuple_t final[1];
        final[0].idx = SA_B[sendcounts[world_rank] - 1].B;
        memcpy(final[0].seq, recvbuf[sendcounts[world_rank] - 1].seq, sizeof(recvbuf[0].seq));
        MPI_Send(&final, 1, MPI_TUPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        //receiving last element from neighbor and checking for overlapping buckets
        tuple_t curr[1];
        MPI_Recv(&curr, 1, MPI_TUPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        int i = 0;
        while (char_array_comp(curr[0].seq, recvbuf[i].seq, K_SIZE))
        {
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    std::vector<tuple_t>().swap(recvbuf);
    MPI_Barrier(MPI_COMM_WORLD);
   
    //################################################//
    // Main for loop to perfrom prefix doubling.      //
    //################################################//
    
    for(int h=K_SIZE;h<=string_length-K_SIZE+1;h*=2)
    {

    //buffer to send out to other processors
    std::vector<tuple_ISA> sendbufISA(sendcounts[world_rank]);

    int displace_ISA[world_size];
    int current_index_ISA=0;
    int counts_filled[world_size];//used to fill the sendbuf ISA
    for(int i=0;i<world_size;i++){
        displace_ISA[i]=current_index_ISA;
        current_index_ISA+=counts_bucket[i];
        counts_filled[i]=0;
    }

    //checking to which processor element has to go 
    for(int j=0;j<sendcounts[world_rank];j++){
        int id= bucket_id(displs, SA_B[j], world_size);
        sendbufISA[displace_ISA[id]+counts_filled[id]]=SA_B[j];
        counts_filled[id]++;
    }

    //Sending out to other processors
    std::vector<tuple_ISA> recvbuf_ISA = probing_alltoallv(sendbufISA, displace_ISA, sendcounts[world_rank], world_size, counts_bucket, MPI_COMM_WORLD, world_rank, MPI_TUPLE_ISA);
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<tuple_ISA>().swap(sendbufISA);

    //declaring B and reordering to string order
    std::vector<int> B(sendcounts[world_rank]);
    reorder_to_stringorder(B,recvbuf_ISA,sendcounts[world_rank],displs[world_rank]);
    std::vector<tuple_ISA>().swap(recvbuf_ISA);

    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displs[i]+sendcounts[i]-1;
    }

    std::vector<int> B2(sendcounts[world_rank]);
    std::copy(B.begin(),B.begin() + sendcounts[world_rank],B2.begin());

    //SHIFTING
    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displs, sendcounts[world_rank], offsets);


    //TRIPLES
    std::vector<triple_t> triple_arr(sendcounts[world_rank]);
    create_triple(B,B2,sendcounts[world_rank],displs[world_rank], triple_arr);
  
    std::vector<int>().swap(B);
    std::vector<int>().swap(B2);

    //SORT LOCALLY
    t_sort(triple_arr, sendcounts[world_rank]);

    MPI_Barrier(MPI_COMM_WORLD);

    //SORT GLOBALLY
    std::vector<triple_t> global_result_triplet;

    if (world_rank == MASTER){
        global_result_triplet = typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
    }
    else
    {
        typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
    }
    std::vector<triple_t>().swap(triple_arr);

    MPI_Barrier(MPI_COMM_WORLD);


    // Set up triple struct array of ints 
    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    MPI_Aint displacements_triple[3] = {0, sizeof(int), 2*sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);

    //SCATTER
    std::vector<triple_t> recvbuf_triplets(sendcounts[world_rank]);
    MPI_Scatterv(&global_result_triplet[0], sendcounts, displs, MPI_TRIPLE_STRUCT, &recvbuf_triplets[0], sendcounts[world_rank], MPI_TRIPLE_STRUCT, 0, MPI_COMM_WORLD);

    std::vector<triple_t>().swap(global_result_triplet);

    //REBUCKET
    //counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    rebucketing(SA_B,recvbuf_triplets, sendcounts[world_rank], displs,world_rank,world_size,counts_bucket);
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank < world_size - 1)
    {
        triple_t final[1];
        final[0].b= recvbuf_triplets[sendcounts[world_rank] - 1].b;
        final[0].b2= recvbuf_triplets[sendcounts[world_rank] - 1].b2;
        final[0].idx= SA_B[sendcounts[world_rank] - 1].B;
        MPI_Send(&final, 1, MPI_TRIPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        triple_t curr[1];
        MPI_Recv(&curr, 1, MPI_TRIPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        int i = 0;
        while (curr[0].b == recvbuf_triplets[i].b and curr[0].b2 == recvbuf_triplets[i].b2)
        {
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }

    std::vector<triple_t>().swap(recvbuf_triplets);

    // check singleton 
    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, sendcounts[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton_global;
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&singleton_global,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

    if(singleton_global){

        MPI_Barrier(MPI_COMM_WORLD);
        //get time of algo
        end = MPI_Wtime();
        MPI_Finalize();

        if (world_rank == 0) {
            //print time elapsed on master node
            std::cout <<  "*" << end-start << "*" << std::endl;
          }

        return 0;
    }

    }
    std::cout << "Singleton Done: " << MPI_Wtime() << "|" << world_rank << std::endl;
    return 0;
    
}
