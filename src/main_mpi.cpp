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
    //setting up timing variables
    double start, end;
    MPI_Init(NULL, NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    // Starting timer
    start = MPI_Wtime();

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


    //std::cout << "Data Reading" << std::endl;


    //std::ifstream in("./input.txt");
    //std::string input((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    //int string_length = input.length();

    //reading in file to get the string length
    std::ifstream is;
    is.open ("./input.txt", std::ios::binary );
    is.seekg (0, std::ios::end);
    int string_length = is.tellg();
    is.close();
    //std::cout<<string_length<<std::endl;

    //std::cout << "Data Read: Done!" << std::endl;
    //char input[]="MISSISSIPPI$";

    int nmin = string_length / world_size;   //min step size to send
    int nextra = string_length % world_size; //remainder if not divisible by nprocessors

    int k = 0;
    int sendnums[world_size]; //size to send out to process i
    int displace[world_size];  //displacement from start
    //maybe do this world_size-1 to not include 0
    for (int i = 0; i < world_size; i++)
    {
        if (i < nextra) sendnums[i] = nmin + 1;
        else sendnums[i] = nmin;
        displace[i] = k;
        k = k + sendnums[i];
    }


    std::fstream File("./input.txt", std::ios::in | std::ios::out );
    File.seekg(displace[world_rank], std::ios::beg);
    char* input;
    int size;
    if(world_rank==world_size-1){
        //input=(char*)malloc(sizeof(char)*sendnums[world_rank]);
        input=new char[sendnums[world_rank]];
        File.read(input, sendnums[world_rank]);
        size=sendnums[world_rank]-K+1;
        std::cout<<"size"<<size<<std::endl;
        //print_char_array(input,sendnums[world_rank]);
    }
    else{
        //input=(char*)malloc(sizeof(char)*(sendnums[world_rank]+K-1));
        input=new char[sendnums[world_rank]+K-1];
        File.read(input, sendnums[world_rank]+K-1);
        size=sendnums[world_rank];
        //print_char_array(input,sendnums[world_rank]+K-1);
    }
    //F[5] = 0;
    //print_char_array(input,sendnums[world_rank]);
    File.close();


    //int string_length = input.length();

    std::vector<tuple_t> global_result_kmers;


  std::cout << "Sorting." << std::endl;

   

    if (world_rank == MASTER)
    {

        //Send out overlapping substrings to workers

        /*
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Send(&input[displace[rank]], sendnums[rank] + K - 1, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        }
        
        */

        //std::cout << "Allocating kmers!" << std::endl;
        //int size = sendnums[0] - K + 1;

        std::vector<tuple_t> kmers(size);
        //tuple_t kmers[size];

        std::cout <<"world_rank"<<world_rank<<"Getting kmers!" << std::endl;

        get_kmers_adapt(input, K, kmers, size,displace[0]);
        //free(input);
        delete[] input;
        std::cout<<"sorting stuff"<<std::endl;
        t_sort(kmers, size);
        //tuple_t_print(kmers, size);

        //std::cout << "Begin global sort!" << std::endl;

        global_result_kmers.resize(string_length - K + 1);
        global_result_kmers = typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
        kmers.clear();
        kmers.shrink_to_fit();
        //t_print(global_result_kmers,string_length - K + 1);

    }
    else
    {
        //int size = sendnums[world_rank]+K-1;

        //char sub_input[size];

        //int recieved_size;
        //MPI_Status status;
        //MPI_Recv(&sub_input[0], size,MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, NULL);
        //MPI_Get_count(&status, MPI_CHAR, &recieved_size);

        //std::cout << "Allocating kmers!" << world_rank << std::endl;
        std::vector<tuple_t> kmers(size);

        //std::cout <<"world_rank"<<world_rank<<"Getting kmers!" << std::endl;
        get_kmers_adapt(input, K, kmers, size,displace[world_rank]);
        //free(input);
        delete[] input;
        //std::cout<<"sorting stuff"<<world_rank<<std::endl;


        // sort the tuples lexicographically
        t_sort(kmers, size);
        //MPI_Barrier(MPI_COMM_WORLD);//should we insert this barrier or not?


         //std::cout << "Begin global sort on worker!"<<world_rank << std::endl;

        //Global sort
        typename_t_sort(log2(world_size), world_rank, kmers, size, MPI_COMM_WORLD);
        kmers.clear();
        kmers.shrink_to_fit();
        //std::cout << "Soting done on " << world_rank<< std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //std::cout<<"fertig"<<std::endl;
    

    //if (world_rank == 0)
     //   //std::cout << (*global_result_kmers).size()<<std::endl;

    //std::cout << "Sorting done!" << std::endl;

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
    //////std::cout<<"before scattering"<<std::endl;
    MPI_Scatterv(&global_result_kmers[0], sendcounts, displs, MPI_TUPLE_STRUCT,
     &recvbuf[0], sendcounts[world_rank], MPI_TUPLE_STRUCT, 0, MPI_COMM_WORLD);

    //tuple_ISA SA_B[sendcounts[world_rank]];//array of tuples on this processor
    std::vector<tuple_ISA> SA_B(sendcounts[world_rank]);

    int counts_bucket[world_size];//counts of tuples to send out to each processor
    //according to their index for later ISA creation
    for(int i=0;i<world_size;i++){
        counts_bucket[i]=0;
    }

    //std::cout << "Rebucketing" << std::endl;

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
        //////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        //////std::cout<<curr[0].idx<<std::endl;
        while (char_array_comp(curr[0].seq, recvbuf[i].seq, K))
        {
            //////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    recvbuf.clear();
    recvbuf.shrink_to_fit();
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

    ////std::cout << "Getting Bucket id";

    //////std::cout << "world rank after sendbuf made"<<world_rank << std::endl;

    for(int j=0;j<sendcounts[world_rank];j++){
        //int id= bucket_id(displs, SA_B[j], world_size);
        int id= bucket_id(displs, SA_B[j], world_size);
        sendbufISA[displace_ISA[id]+counts_filled[id]]=SA_B[j];
        counts_filled[id]++;
    }

    std::cout << "Pro ing all to allv" <<std::endl;





    //this array must be updated after every for loop
    //std::vector<tuple_ISA> recvbuf_ISA(sendcounts[world_rank]);
    std::vector<tuple_ISA> recvbuf_ISA = probing_alltoallv(sendbufISA, displace_ISA, sendcounts[world_rank], world_size, counts_bucket, MPI_COMM_WORLD, world_rank, MPI_TUPLE_ISA);
    //problem: displacements in recvbuf not known. Could we probe?
    //MPI_Alltoallv(&sendbufISA,counts_bucket,displace_ISA,MPI_TUPLE_ISA,recvbuf_ISA,)
    //sendbufISA.clear();
    //sendbufISA.shrink_to_fit();
    MPI_Barrier(MPI_COMM_WORLD);
    sendbufISA.clear();
    sendbufISA.shrink_to_fit();

    std::cout << "Reordering" << std::endl;

    std::vector<int> B(sendcounts[world_rank]);
    reorder_to_stringorder(B,recvbuf_ISA,sendcounts[world_rank],displs[world_rank]);
    recvbuf_ISA.clear();
    recvbuf_ISA.shrink_to_fit();

    //print_int_array(B,sendcounts[world_rank]);

    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displs[i]+sendcounts[i]-1;
    }

    std::vector<int> B2(sendcounts[world_rank]);
    std::copy(B.begin(),B.begin() + sendcounts[world_rank],B2.begin());

    if(world_rank==0){

    std::cout << "Reordered!" << std::endl;


    std::cout << "Shifting" << std::endl;
    }

    //SHIFTING
    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displs, sendcounts[world_rank], offsets);

    std::cout <<world_rank<<"Shifting done!" << std::endl;
    //TRIPLES

    //print_int_array(B2,sendcounts[world_rank]);
    

    std::vector<triple_t> triple_arr(sendcounts[world_rank]);
    create_triple(B,B2,sendcounts[world_rank],displs[world_rank], triple_arr);
    B.clear();
    B.shrink_to_fit();
    B2.clear();
    B2.shrink_to_fit();

    //////std::cout << "-----" << std::endl;


    //std::cout << "Sorting Triples" << std::endl;

    //SORT LOCALLY
    t_sort(triple_arr, sendcounts[world_rank]);

    MPI_Barrier(MPI_COMM_WORLD);

    //SORT GLOBALLY
    std::vector<triple_t> global_result_triplet;

    if (world_rank == MASTER){
        global_result_triplet = typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
        //t_print(global_result_triplet,string_length-K+1);
    }
    else
    {
        typename_t_sort(log2(world_size), world_rank, triple_arr, sendcounts[world_rank], MPI_COMM_WORLD);
    }
    triple_arr.clear();
    triple_arr.shrink_to_fit();
    MPI_Barrier(MPI_COMM_WORLD);

    

    //std::cout << "Triples sorted!" << std::endl;

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
        //////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K)<<std::endl;
        //////std::cout<<curr[0].idx<<std::endl;
        while (curr[0].b == recvbuf_triplets[i].b and curr[0].b2 == recvbuf_triplets[i].b2)
        {
            //////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx; //CORRECT i was not there
            i++;
        }
    }
    recvbuf_triplets.clear();
    recvbuf_triplets.shrink_to_fit();

    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, sendcounts[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton_global;
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&singleton_global,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

    if(singleton_global){
        std::cout<<"singleton is global"<<std::endl;
        std::vector<tuple_ISA> final_SA(string_length-K+1);
        MPI_Gatherv(&SA_B[0],sendcounts[world_rank],MPI_TUPLE_ISA,&final_SA[0],sendcounts,displs,MPI_TUPLE_ISA,0,MPI_COMM_WORLD);
        std::cout<<"just gathered lol"<<std::endl;
        if(world_rank==0){
            debug_tuple_print(final_SA,string_length-K+1);
            std::cout<<"please print this"<<std::endl;

        }

        MPI_Barrier(MPI_COMM_WORLD);
        //get time of algo

        end = MPI_Wtime();
        std::cout<<"Am I done??"<<std::endl;

        MPI_Finalize();
        std::cout<<"Am I done??"<<std::endl;
        if (world_rank == 0) {
             //print time elapsed on master node

            std::cout <<  "*" << end-start << "*" << std::endl;
          }

        return 0;

        //exit(1);
    }

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    // Finalize the MPI environment.
    
    //MPI_Finalize();
    return 0;
    
}
