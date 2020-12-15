#include <iostream>
#include <string>
#include <mpi.h>
#include <math.h>
#include "helper.cpp"
#include <fstream>
#include <vector>

#define MASTER 0

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    //SETUP TUPLE STRUCT
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K};
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

    
    /*
    char* input;
    if (argc > 1){
        input = argv[1];
        //K = argv[2];
    }
    */
    

    //std::cout << "Data Reading" << std::endl;

    
    std::ifstream in("./input.txt");
    std::string input((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    int string_length = input.length();

    std::cout << "Data Read: Done!" << std::endl;

    //char input[] = "agaagccagtactgcgacaaaggtaggacatggcgttgcaccaaatcagtaccggctccacaataattacaccatagggcaccgctatccgcgtgcgtca$";
    //char input[] = "ABCDEFGHAHAHABCDEFGHAHAHABCDEFGHAHAH$";
    //char input[] = "AAAASBBBL"; //overlapping string of k=2 and p=4
    //char input[]="NJSDJSNMCXNCJDKHFUIAAAJ";
    //char input[]="RAYASTOYANOVA";
    //char input[]="MISSISSIPPI$";

    int nmin = string_length / world_size;   //min step size to send
    int nextra = string_length % world_size; //remainder if not divisible by nprocessors

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

    std::vector<tuple_t> global_result_kmers;


   std::cout << "Sorting." << std::endl;

    if (world_rank == MASTER)
    {

        //Send out overlapping substrings to workers

        
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Send(&input[displace[rank]], sendnums[rank] + K - 1, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        }
        
        std::cout << "Allocating kmers!" << std::endl;
        int size = sendnums[0] - K + 1;
        
        std::vector<tuple_t> kmers(size);
        //tuple_t kmers[size];

        std::cout << "Getting kmers!" << std::endl;

        get_kmers_adapt(&input[displace[0]], K, kmers, size,displace[0]);
        t_sort(kmers, size);
        //tuple_t_print(kmers, size);

        std::cout << "Begin global sort!" << std::endl;

        global_result_kmers.resize(string_length - K + 1);
        global_result_kmers = typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
        
    }
    else
    {
        int size = sendnums[world_rank]+K-1;
 
        char sub_input[size];

        //int recieved_size;
        //MPI_Status status;
        MPI_Recv(&sub_input[0], size,MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, NULL);
        //MPI_Get_count(&status, MPI_CHAR, &recieved_size);

        std::cout << "Allocating kmers!" << world_rank << std::endl;
        std::vector<tuple_t> kmers(size - K + 1);

        std::cout << "Getting kmers!" << std::endl;
        get_kmers_adapt(sub_input, K, kmers, size - K + 1,displace[world_rank]);
   
        
        // sort the tuples lexicographically
        t_sort(kmers, size - K + 1);
        //MPI_Barrier(MPI_COMM_WORLD);//should we insert this barrier or not?
    

         std::cout << "Begin global sort on worker!"<<world_rank << std::endl;

        //Global sort
        typename_t_sort(log2(world_size), world_rank, kmers, size - K + 1, MPI_COMM_WORLD);
        std::cout << "Soting done on " << world_rank<< std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    //if (world_rank == 0)
     //   std::cout << (*global_result_kmers).size()<<std::endl;

    std::cout << "Sorting done!" << std::endl;

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


    std::vector<tuple_t> recvbuf(sendcounts[world_rank]);
    ////std::cout<<"before scattering"<<std::endl;
    MPI_Scatterv(&global_result_kmers[0], sendcounts, displs, MPI_TUPLE_STRUCT,
     &recvbuf[0], sendcounts[world_rank], MPI_TUPLE_STRUCT, 0, MPI_COMM_WORLD);
    
    //tuple_ISA SA_B[sendcounts[world_rank]];//array of tuples on this processor
    std::vector<tuple_ISA> SA_B(sendcounts[world_rank]);

    int counts_bucket[world_size];//counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    std::cout << "Rebucketing" << std::endl;
    
    rebucketing(SA_B, recvbuf, sendcounts[world_rank], displs,counts_bucket,world_rank, world_size);
    

    
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
        ////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        ////std::cout<<curr[0].idx<<std::endl;
        while (char_array_comp(curr[0].seq, recvbuf[i].seq, K))
        {
            ////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    std::cout << "Rebucketing done!" << std::endl;


    //t_print(SA_B, sendcounts[world_rank]); Correct

    //########################################################
    //########################################################
    //########################################################
    //########################################################


    for(int h=K;h<=string_length-K+1;h*=2)
    {
    
    std::vector<tuple_ISA> sendbufISA(sendcounts[world_rank]);

    int displace_ISA[world_size];
    int current_index_ISA=0;
    int counts_filled[world_size];//used to fill the sendbuf ISA
    for(int i=0;i<world_size;i++){
        displace_ISA[i]=current_index_ISA;
        current_index_ISA+=counts_bucket[i];
        counts_filled[i]=0;
    }

    //std::cout << "Getting Bucket id";

    ////std::cout << "world rank after sendbuf made"<<world_rank << std::endl;

    for(int j=0;j<sendcounts[world_rank];j++){
        //int id= bucket_id(displs, SA_B[j], world_size);
        int id= bucket_id(displs, SA_B[j], world_size);
        sendbufISA[displace_ISA[id]+counts_filled[id]]=SA_B[j];
        counts_filled[id]++;    
    }

    //std::cout << "Pro ing all to allv" <<std::endl;


  


    //this array must be updated after every for loop
    std::vector<tuple_ISA> recvbuf_ISA = probing_alltoallv(sendbufISA, displace_ISA, sendcounts[world_rank], world_size, counts_bucket, MPI_COMM_WORLD, world_rank, MPI_TUPLE_ISA);
    

    
    MPI_Barrier(MPI_COMM_WORLD);
    
    std::cout << "Reordering" << std::endl;
    
    std::vector<int> B(sendcounts[world_rank]);
    reorder_to_stringorder(B,recvbuf_ISA,sendcounts[world_rank],displs[world_rank]);

    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displs[i]+sendcounts[i]-1;
    }

    std::vector<int> B2(sendcounts[world_rank]);
    std::copy(B.begin(),B.begin() + sendcounts[world_rank],B2.begin());

    std::cout << "Reordered!" << std::endl;


    std::cout << "Shifting" << std::endl;

    //SHIFTING 
    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displs, sendcounts[world_rank], offsets);

    std::cout << "Shifting done!" << std::endl;
    //TRIPLES 

    std::vector<triple_t> triple_arr(sendcounts[world_rank]);
    create_triple(B,B2,sendcounts[world_rank],displs[world_rank], triple_arr);
    ////std::cout << "-----" << std::endl;


    std::cout << "Sorting Triples" << std::endl;
    
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
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Triples sorted!" << std::endl;

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
        ////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        ////std::cout<<curr[0].idx<<std::endl;
        while (curr[0].b == recvbuf_triplets[i].b and curr[0].b2 == recvbuf_triplets[i].b2)
        {
            ////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx; //CORRECT i was not there
            i++;
        }
    }


    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, sendcounts[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton_global;
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&singleton_global,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

    if(singleton_global){
        std::vector<tuple_ISA> final_SA(string_length-K+1);
        MPI_Gatherv(&SA_B[0],sendcounts[world_rank],MPI_TUPLE_ISA,&final_SA[0],sendcounts,displs,MPI_TUPLE_ISA,0,MPI_COMM_WORLD);
        if(world_rank==0){
            //int finalarr[]={11,10,7,4,1,0,9,8,6,3,5,2};

            debug_tuple_print(final_SA,string_length-K+1);

        }
        MPI_Finalize();
        return 0;

        //exit(1);
    }

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    // Finalize the MPI environment.
    //MPI_Finalize();
    return 0;
}