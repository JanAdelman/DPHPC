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
    if(argv[1]) INPUT_PATH = argv[1];
    std::cout << INPUT_PATH << std::endl;	
    //showMemUsage("after initialization", world_rank);

    //SETUP TUPLE STRUCT
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, 1};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
    MPI_Datatype types_ISA[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths_ISA, displacements_ISA, types_ISA, &MPI_TUPLE_ISA);
    MPI_Type_commit(&MPI_TUPLE_ISA);

    //showMemUsage("after initialization of datatypes", world_rank);

    

    //MAIN LOGIC


    //std::cout << "Data Reading: " << MPI_Wtime() << "|" << world_rank << std::endl;


    //std::ifstream in("./input.txt");
    //std::string input((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    //int string_length = input.length();

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
    //maybe do this world_size-1 to not include 0
    for (int i = 0; i < world_size; i++)
    {
        if (i < kextra)
            sendk[i] = kmin + 1;
        else
            sendk[i] = kmin;
        displsk[i] = now;
        now = now + sendk[i];
    }

    /*
    for(int i=0;i<world_size;i++){
        std::cout<<"sendk: "<<sendk[i]<<std::endl;
        std::cout<<"displsk:"<<displsk[i]<<std::endl;
    }
    */

    int sendnums[world_size]; //size to send out to process i
    int displace[world_size];  //displacement from start
    //maybe do this world_size-1 to not include 0
    sendnums[0]=sendk[0]+K_SIZE-1;
    displace[0]=0;
    for (int i = 1; i < world_size; i++)
    {
        sendnums[i]=sendk[i]+K_SIZE-1;
        displace[i]=displace[i-1]+sendnums[i-1]-K_SIZE+1;
    }

    //for(int i=0;i<world_rank;i++){
    //    std::cout<<sendnums[i]<<std::endl;
    //}

    //showMemUsage("getting displ and sendnums", world_rank);


    std::fstream File(INPUT_PATH, std::ios::in | std::ios::out );
    File.seekg(displace[world_rank], std::ios::beg);
    char* input=new char[sendnums[world_rank]];
    File.read(input, sendnums[world_rank]);
    int size=sendk[world_rank];
    File.close();
    //F[5] = 0;
    //print_char_array(input,sendnums[world_rank]);
    //std::cout<<world_rank<<"world rank "<<size<<"size"<<std::endl;
    
    //showMemUsage("After opening file", world_rank);

    //std::cout << "Data Read: " << MPI_Wtime() << "|" << world_rank << std::endl;


    //std::cout<<"wr: "<<world_rank<<"size: "<<size<<std::endl;
    std::vector<tuple_t<unsigned long long int, int>> kmers(size);
    //std::cout<<world_rank<<"Made vector"<<std::endl;
    get_kmers_adapt(input, K_SIZE, kmers, size,displace[world_rank]);
    MPI_Barrier(MPI_COMM_WORLD);
    //showMemUsage("After getting kmers", world_rank);
    delete[] input;
    //showMemUsage("After deleting input", world_rank);
    //std::cout<<"this is the capacity before lol:"<<world_rank<<"wr was this"<<kmers.capacity()<<std::endl;
    samplesort<unsigned long long int,int>(kmers,MPI_COMM_WORLD);
    //std::cout<<"this is the capacity after lol:"<<world_rank<<"wr was this"<<kmers.capacity()<<std::endl;
    //showMemUsage("samplesort", world_rank);
    //tup_t_print<unsigned long long int, int>(kmers,size,world_rank);


    MPI_Barrier(MPI_COMM_WORLD);


    std::vector<tuple_ISA<int>> SA_B(size);

    rebucketing_encode(SA_B, kmers, size, displsk,world_rank, world_size);


    MPI_Barrier(MPI_COMM_WORLD);
    
    

    if (world_rank < world_size - 1)
    {
        tuple_t<unsigned long long int,int> final[1];
        final[0].idx = SA_B[size - 1].B;
        final[0].seq=kmers[size - 1].seq;
        //memcpy(final[0].seq, recvbuf[sendcounts[world_rank] - 1].seq, sizeof(recvbuf[0].seq));
        //final[0].seq=recvbuf[sendcounts[world_rank]-1].seq;
        MPI_Send(&final, 1, MPI_TUPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        tuple_t<unsigned long long int, int> curr[1];
        MPI_Recv(&curr, 1, MPI_TUPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        //tuple_t_print(curr, K_SIZE);
        int i = 0;
        ////////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K_SIZE)<<std::endl;
        ////////std::cout<<curr[0].idx<<std::endl;
        while (curr[0].seq==kmers[i].seq)
        {
            ////////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx;
            i++;
        }
    }
    std::vector<tuple_t<unsigned long long int, int>>().swap(kmers);
    MPI_Barrier(MPI_COMM_WORLD);

    


    //std::cout << "Rebucketing done: " << MPI_Wtime() << "|" << world_rank << std::endl;


    int offsets[world_size];
    for(int i=0;i<world_size;i++){
        offsets[i]=displsk[i]+sendk[i]-1;
    }

    for(int h=K_SIZE;h<=string_length-K_SIZE+1;h*=2)
    {

    
    reorder_string<int>(SA_B,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<int> B(size);
    make_B(SA_B,B,size);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //std::cout << "got offsets " << MPI_Wtime() << "|" << world_rank << std::endl;

    
    std::vector<int> B2(B);

    //std::copy(B.begin(),B.begin() + sendcounts[world_rank],B2.begin());
    //std::cout << "copied stuff " << MPI_Wtime() << "|" << world_rank << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    //std::cout << "Reorder to String Order Done: " << MPI_Wtime() << "|" << world_rank << std::endl;

    //std::cout << "Shifting: " << MPI_Wtime() << "|" << world_rank << std::endl;
    //SHIFTING

    naive_shift(B2, h, MPI_COMM_WORLD, world_rank, world_size, displsk, size, offsets);

    //std::cout << "Shifting done: " << MPI_Wtime() << "|" << world_rank << std::endl;

    //TRIPLES
    //std::cout<<"h is now"<<h<<"on world rank:"<<world_rank<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    

    std::vector<triple_t<int>> triple_arr(size);
    create_triple(B,B2,size,displsk[world_rank], triple_arr);
    std::vector<int>().swap(B);
    std::vector<int>().swap(B2);


    //std::cout << "Sort Triples: " << MPI_Wtime() << "|" << world_rank << std::endl;



    //SORTing triples
    triple_sort<int>(triple_arr,MPI_COMM_WORLD);
    

    MPI_Barrier(MPI_COMM_WORLD);

    

    //std::cout << "Sort Triples Done: " << MPI_Wtime() << "|" << world_rank << std::endl;
    

    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    MPI_Aint displacements_triple[3] = {0, sizeof(int), 2*sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);

    //SCATTER

    //REBUCKET


    //std::cout << "Rebucketing: " << MPI_Wtime() << "|" << world_rank << std::endl;


    //tuple_ISA SA_B[sendcounts[world_rank]];
    rebucketing(SA_B,triple_arr, size, displsk,world_rank,world_size);

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank < world_size - 1)
    {
        triple_t<int> final[1];
        final[0].b= triple_arr[size - 1].b;
        final[0].b2= triple_arr[size - 1].b2;
        final[0].idx= SA_B[size - 1].B;
        //final[0].seq=recvbuf[sendcounts[world_rank]-1].seq;
        MPI_Send(&final, 1, MPI_TRIPLE_STRUCT, world_rank + 1, 0, MPI_COMM_WORLD);
    }
    if (world_rank != 0)
    {
        triple_t<int> curr[1];
        MPI_Recv(&curr, 1, MPI_TRIPLE_STRUCT, world_rank - 1, 0, MPI_COMM_WORLD, NULL);
        //tuple_t_print(curr, K_SIZE);
        int i = 0;
        ////////std::cout<<char_array_comp(curr[0].seq,recvbuf[i].seq,K_SIZE)<<std::endl;
        ////////std::cout<<curr[0].idx<<std::endl;
        while (curr[0].b == triple_arr[i].b and curr[0].b2 == triple_arr[i].b2)
        {
            ////////std::cout<<"while_loop"<<std::endl;
            SA_B[i].B = curr[0].idx; //CORRECT i was not there
            i++;
        }
    }
     std::vector<triple_t<int>>().swap(triple_arr);


    bool singleton = all_singleton(SA_B,MPI_COMM_WORLD, world_rank, world_size, size);
    MPI_Barrier(MPI_COMM_WORLD);

    bool singleton_global;
    MPI_Reduce(&singleton, &singleton_global, 1, MPI_C_BOOL, MPI_LAND, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&singleton_global,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    //std::cout<<"h: "<<h<<"all singleton"<<singleton_global<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if(singleton_global){
        
	std::vector<tuple_ISA<int>> final_SA;
        if (world_rank == 0)
            final_SA.resize(string_length-K_SIZE+1);
        
	MPI_Gatherv(&SA_B[0],size,MPI_TUPLE_ISA,&final_SA.front(),sendk,displsk,MPI_TUPLE_ISA,0,MPI_COMM_WORLD);

        if(world_rank==0){
            debug_tuple_print(final_SA,string_length-K_SIZE+1);
        }
	
        MPI_Barrier(MPI_COMM_WORLD);
        //get time of algo
	
	mytime = MPI_Wtime() - mytime;
        end = MPI_Wtime();

	double maxtime, mintime, avgtime;
	MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
	MPI_Reduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        
	MPI_Finalize();

        if (world_rank == 0) {
             //print time elapsed on master node
	    avgtime /= world_size;
	    	
            std::cout << "Time of the algo:  "<< "*" << end-start << "*" << std::endl;
	    std::cout << "min_time: " << "*" <<  mintime << "*" << std::endl;
	    std::cout << "max_time: " << "*" <<  maxtime << "*" <<std::endl;
            std::cout << "avergage_time: "<< "*" << avgtime << "*" << std::endl;
          }

        return 0;

        //exit(1);
    }

    }
    //std::cout << "Singleton Done: " << MPI_Wtime() << "|" << world_rank << std::endl;


    //MPI_Barrier(MPI_COMM_WORLD);
    // Finalize the MPI environment.
    
}
