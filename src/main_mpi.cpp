#include <iostream>
#include <string>
#include <mpi.h>
#include <math.h>
#include "helper.cpp"
#include <fstream>
#include <vector>
#include "mxx/include/mxx/comm.hpp"
#include "mxx/include/mxx/sort.hpp"
#include "mxx/include/mxx/datatypes.hpp"

#define MASTER 0

std::string INPUT_PATH = "./data/input.txt";

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    // setting up timing variables 
    double start, end, mytime;

    //Start timing 
    start = MPI_Wtime();
    mytime = MPI_Wtime(); //timing for min/max/average 

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
    // Command line argument for the path to the data 
    if(argv[1]) INPUT_PATH = argv[1];	
    
    //SETUP TUPLE STRUCT (Costum MPI Datatype)
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, 1};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);
    
    // Set up of Tuple_ISA
    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
    MPI_Datatype types_ISA[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths_ISA, displacements_ISA, types_ISA, &MPI_TUPLE_ISA);
    MPI_Type_commit(&MPI_TUPLE_ISA);

    //reading in file to get the string length   
    std::ifstream is;
    is.open (INPUT_PATH, std::ios::binary );
    is.seekg (0, std::ios::end);
    int string_length = is.tellg();
    is.close();
    
    int kmin = (string_length - K_SIZE + 1) / world_size;   //min step size to send
    int kextra = (string_length - K_SIZE + 1) % world_size; //remainder if not divisible by nprocessors
    int now = 0;
    int sendk[world_size]; //size to send out to process i
    int displsk[world_size];     //displacement from start

    for (int i = 0; i < world_size; i++)
    {
        if (i < kextra)
            sendk[i] = kmin + 1;
        else
            sendk[i] = kmin;
        displsk[i] = now;
        now = now + sendk[i];
    }


    int sendnums[world_size]; //size to send out to process i
    int displace[world_size];  //displacement from start
    sendnums[0]=sendk[0]+K_SIZE-1;
    displace[0]=0;
    for (int i = 1; i < world_size; i++)
    {
        sendnums[i]=sendk[i]+K_SIZE-1;
        displace[i]=displace[i-1]+sendnums[i-1]-K_SIZE+1;
    }


    std::fstream File(INPUT_PATH, std::ios::in | std::ios::out );
    File.seekg(displace[world_rank], std::ios::beg);
    char* input=new char[sendnums[world_rank]];
    File.read(input, sendnums[world_rank]);
    int size=sendk[world_rank];
    File.close();
   
    // set up kmers vector 
    std::vector<tuple_t<unsigned long long int, int>> kmers(size);
    
    // get kmers on each process and the string is also encoded 
    get_kmers_adapt(input, K_SIZE, kmers, size,displace[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);
	
    // delete input string 
    delete[] input;
   
    // sort the kmers 
    samplesort<unsigned long long int,int>(kmers,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // create vector SA_B 
    std::vector<tuple_ISA<int>> SA_B(size);
	
    // rebucket on each process
    rebucketing_encode(SA_B, kmers, size, displsk,world_rank, world_size);
    MPI_Barrier(MPI_COMM_WORLD);
        
    // get portion to send expect (last process only receives) 
    if (world_rank < world_size - 1)
    {
        tuple_t<unsigned long long int,int> final[1];
        final[0].idx = SA_B[size - 1].B;
        final[0].seq=kmers[size - 1].seq;
        MPI_Send(&final, 1, MPI_TUPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    // receive the data if you are not process 0 
    if (world_rank != 0)
    {
        tuple_t<unsigned long long int, int> curr[1];
        MPI_Recv(&curr, 1, MPI_TUPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        int i = 0;
        while (curr[0].seq==kmers[i].seq)
        {
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    // deallocate kmers vector 
    std::vector<tuple_t<unsigned long long int, int>>().swap(kmers);
    MPI_Barrier(MPI_COMM_WORLD);

    // get offsets 
    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displsk[i]+sendk[i]-1;
    }

    //##########################################//
    // Main Prefix doubling loop starts here.   //
    //##########################################//
	
    for(int h=K_SIZE;h<=string_length-K_SIZE+1;h*=2)
    {

    // perform reorder to string order on SA_B
    reorder_string<int>(SA_B,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // set up vector B 
    std::vector<int> B(size);
    make_B(SA_B,B,size);
    MPI_Barrier(MPI_COMM_WORLD);
 
    // set up vector B2
    std::vector<int> B2(B);  
    MPI_Barrier(MPI_COMM_WORLD);


    //SHIFTING by h
    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displsk, size, offsets);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // set up triple vector 
    std::vector<triple_t<int>> triple_arr(size);
    create_triple(B,B2,size,displsk[world_rank], triple_arr);
    // dellocate vectors B and B2	     
    std::vector<int>().swap(B);
    std::vector<int>().swap(B2);

    //Sorting triples
    triple_sort<int>(triple_arr,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // set up costum MPI datastructre for the vectors of triples 
    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    MPI_Aint displacements_triple[3] = {0, sizeof(int), 2*sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);

    //SCATTER
    // rebucket 
    rebucketing(SA_B,triple_arr, size, displsk,world_rank,world_size);
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank < world_size - 1)
    {
        triple_t<int> final[1];
        final[0].b= triple_arr[size - 1].b;
        final[0].b2= triple_arr[size - 1].b2;
        final[0].idx= SA_B[size - 1].B;
        MPI_Send(&final, 1, MPI_TRIPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        triple_t<int> curr[1];
        MPI_Recv(&curr, 1, MPI_TRIPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        int i = 0;
        while (curr[0].b == triple_arr[i].b and curr[0].b2 == triple_arr[i].b2)
        {
            
            SA_B[i].B = curr[0].idx; 
            i++;
        }
    }
    // remove triple vectoor 	     
    std::vector<triple_t<int>>().swap(triple_arr);

    // check if singleton 
    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, size);
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton_global;
    // checks on each processor if it is singleton     
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&singleton_global,1,MPI_C_BOOL,0,MPI_COMM_WORLD)
    MPI_Barrier(MPI_COMM_WORLD);
	    
    // check if all groups are singleton 
    if(singleton_global){
 	
        MPI_Barrier(MPI_COMM_WORLD);
        //get time of algo
	mytime = MPI_Wtime() - mytime;
        end = MPI_Wtime();
	    
	// compute the min/ max and average time via MPI_Reduce()
	double maxtime, mintime, avgtime;
	MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
	MPI_Reduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        
	MPI_Finalize();

        if (world_rank == 0) {
            //print time elapsed on master node
	    avgtime /= world_size;

            // print the resulting time 
            std::cout << end-start << "," << mintime << "," << maxtime << "," << avgtime << std::endl;
	    

            return 0;

    }

    }

}
