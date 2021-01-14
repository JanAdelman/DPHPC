#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>
#include <type_traits>
#include <ostream>
#include <fstream>
#include "mxx/include/mxx/comm.hpp"
#include "mxx/include/mxx/sort.hpp"
#include "mxx/include/mxx/datatypes.hpp"
#include <bitset>

#define K_SIZE 32

// set up struct of tuples 
template <typename T, typename index_t>
struct tuple_t
{
    T seq;
    index_t idx;
};

// set up triple struct 
template<typename T>
struct triple_t
{
    T idx;
    T b;
    T b2;
};
// set up of inverted suffix array struct 
template <typename T>
struct tuple_ISA{
    T SA;
    T B;
};

namespace mxx {
template <typename T, typename index_t>
MXX_CUSTOM_TEMPLATE_STRUCT(MXX_WRAP_TEMPLATE(tuple_t<T,index_t>), seq, idx);
template <typename ISA>
MXX_CUSTOM_TEMPLATE_STRUCT(MXX_WRAP_TEMPLATE(tuple_ISA<ISA>),SA,B);
template <typename triple>
MXX_CUSTOM_TEMPLATE_STRUCT(MXX_WRAP_TEMPLATE(triple_t<triple>),idx,b,b2);
} 


// get buecket id 
int bucket_id(int* displ, tuple_ISA<int> &SA_B, int world_size){
    for(int i=0;i<world_size-1;i++){
        if (displ[i]<=SA_B.SA && displ[i+1]>SA_B.SA){
            return i;
        }
    }
    return world_size-1;
}

int bucket_id_shift(int* displ, int B, int world_size, int* sendcounts){
    if(B>displ[world_size-1]+sendcounts[world_size-1]-1){
        return -1;
    }
    for(int i=0;i<world_size-1;i++){
        if (displ[i]<=B && displ[i+1]>B){
            return i;
        }
    }
    return world_size-1;
}

// rebucket with the encoded input 
void rebucketing_encode(std::vector<tuple_ISA<int>> &SA_B, std::vector<tuple_t<unsigned long long int,int>> &kmers, size_t size, int* displ, int world_rank, int world_size){
    SA_B[0].B=displ[world_rank];
    SA_B[0].SA=kmers[0].idx;
    for(int i=1;i<size;i++){
        if(kmers[i-1].seq==kmers[i].seq){
            SA_B[i].B=SA_B[i-1].B;
            SA_B[i].SA=kmers[i].idx;
        }
        else{
            SA_B[i].B=displ[world_rank]+i;
            SA_B[i].SA=kmers[i].idx;
        }
    }
}

//rebucketing not encoded 
void rebucketing(std::vector<tuple_ISA<int>> &SA_B,std::vector<triple_t<int>> &input, size_t size, int* displ, int world_rank, int world_size){
    SA_B[0].B=displ[world_rank];
    SA_B[0].SA=input[0].idx;
    for(int i=1;i<size;i++){
        if((input[i-1].b == input[i].b and input[i-1].b2 == input[i].b2)){
            SA_B[i].B=SA_B[i-1].B;
            SA_B[i].SA=input[i].idx;
        }
        else{
            SA_B[i].B=displ[world_rank]+i;
            SA_B[i].SA=input[i].idx;
        }
    }
}


// helper function to print char array 
void print_char_array(const char *input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout <<"   "<< input[i];
    std::cout << std::endl;
}

// helper function to print int array 
void print_int_array(std::vector<int> &input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout << input[i]<<",";
    std::cout << std::endl;
}

// function to encode string of {A,C,T,G} to two bits 
unsigned long long int encode (char* input, int size){
  std::bitset<K_SIZE*2> p;

    static std::map<char, std::tuple<int, int>> table = {
        { 'A', std::make_tuple(0,0) },
        { 'C', std::make_tuple(0,1) },
        { 'G', std::make_tuple(1,0) },
        { 'T', std::make_tuple(1,1) }
    };

    int bit = 0;
    for (int i = size - 1; i >= 0 ; i--){
        p[bit + 1] = std::get<0>(table[input[i]]);
        p[bit] = std::get<1>(table[input[i]]);
        bit += 2;
    }

    return p.to_ullong();
}

// functions to get kmers + encodes the input array to two bits 
void get_kmers_adapt(char *input, const int k, std::vector<tuple_t<unsigned long long int, int>> &kmers, size_t size, int displacement)
{

    for (int i = 0; i < size; i++)
    {
        kmers[i].seq=encode(input+i,k);
        kmers[i].idx = displacement+ i;
    }
}

// compare chars of two arrays, boolean output 
bool char_array_comp(const char* a, const char* b, int size){
    bool same = true;
    for (int i = 0; i < size; i++) {
        if (*(a + i) != *(b + i))
            same = false;
    }
    return same; 
}

// helper function to print tuple structs 
template <typename T, typename index_t>
void tup_t_print(std::vector<tuple_t<T,index_t>> input, size_t size, int world_rank){
    for(int i=0;i<size;i++){
        std::cout<<"seq "<<input[i].seq<<std::endl;
        std::cout<<" idx"<<input[i].idx<<std::endl;
    }
    std::cout<<"world_rank:"<<world_rank<<std::endl;
}
// helper function to print inverted suffix array 
void isa_print(std::vector<tuple_ISA<int>> input, size_t size, int world_rank){
    for(int i=0;i<size;i++){
        std::cout<<"SA "<<input[i].SA<<std::endl;
        std::cout<<" B"<<input[i].B<<std::endl;
    }
    std::cout<<"world_rank:"<<world_rank<<std::endl;
}

// helper function to print triple struct 
void t_print(std::vector<triple_t<int>> &input, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        std::cout << "    B: " << input[i].b << std::endl;
        std::cout << "    B2: " << input[i].b2 << std::endl;
        std::cout << "    SA: " << input[i].idx << std::endl;
    }
    std::cout << "---" << std::endl;
}

// helper function to print triple function
void t_print_flat(const triple_t<int> *input, size_t size, int rank, int target)
{
    std::cout <<rank<<"(from):(to)"<<target<< "->";
    for (int i = 0; i < size; i++)
    {
        std::cout << input[i].b <<","<< input[i].b2<<","<<input[i].idx<<";";//  << std::endl;
    }
    std::cout << std::endl;
}


// herlper function to write output to a file calles result.txt 
void debug_tuple_print(const std::vector<tuple_ISA<int>> input, size_t size)
{
    
    std::ofstream outfile ("./result.txt");
    
    for (int i = 0; i < size; i++)
    {
        outfile << input[i].SA << ",";
    }
    outfile.close();
}

// create B vector 
void make_B(std::vector<tuple_ISA<int>> &input,std::vector<int> &B,size_t size){
    for(int i=0;i<size;i++){
        B[i]=input[i].B;
    }
}

// sample sort function from the mxx library 
template <typename T, typename index_t>
void samplesort(std::vector<tuple_t<T,index_t>>& vec,MPI_Comm comm){
    using TT = tuple_t<T, index_t>;
    auto cmpidx = [](const TT &a, const TT &b){
        if (a.seq<b.seq)
            return true;
        else if (a.seq==b.seq)
            return a.idx<b.idx;
        else 
            return false;
     };
     mxx::sort(vec.begin(),vec.end(),cmpidx,MPI_COMM_WORLD);
}
// rearder to stringorder by sorting 
template <typename T>
void reorder_string(std::vector<tuple_ISA<int>> &input,MPI_Comm comm){
    using TT=tuple_ISA<int>;
    auto cmpidx=[](const TT &a, const TT &b){
        if(a.SA<b.SA){
            return true;
        }
        else{
            return false;
        }
    };
    mxx::sort(input.begin(),input.end(),cmpidx,comm);
}
// sort by uisng mxx:sort and costum compare function 
template <typename T>
void triple_sort(std::vector<triple_t<int>> &input,MPI_Comm comm){
    using TT=triple_t<int>;
    auto cmpidx=[](const TT &a, const TT &b){
        return a.b < b.b || (a.b == b.b && a.b2 < b.b2) || (a.b == b.b && a.b2 == b.b2 && a.idx < b.idx);
    };
    mxx::sort(input.begin(),input.end(),cmpidx,comm);
}

// function to create triple array 
void create_triple(std::vector<int> &B, std::vector<int> &B2, int size,int displ, std::vector<triple_t<int>> &triple_arr){
    for(int i=0;i<size;i++){
        triple_arr[i].b=B[i];
        triple_arr[i].b2=B2[i];
        triple_arr[i].idx=i+displ;  
    }
    
}

// function to perform shifting, shifting over more than one process possible 
void naive_shift(std::vector<int> &input, const int h, MPI_Comm comm, const int world_rank, const int world_size, int *offsets_start, int local_length, int *offsets_end){

            if (world_rank == 0)
            {
                if (h < local_length){
                    std::copy(&input[h], &input[local_length], &input[0]);
                }
            }
            else{
                int shift_start = offsets_start[world_rank] - h;
                int shift_end = offsets_end[world_rank] - h;

                for (int rank = world_rank; rank >= 0; rank--){//Look at previous buckets

                    if ((shift_start >= offsets_start[rank]) and (shift_start <= offsets_end[rank])){//Start is within the offsets of this rank
      
                        if (shift_end >= offsets_start[world_rank]){ //End is also within bounds -> Only send to one porcessor

                            MPI_Send(&input[0], h, MPI_INT, rank, shift_start, comm);
                            std::copy(&input[h], &input[local_length], &input[0]);
                            break;                           
                        }
                        else if (shift_end <= offsets_end[rank])
                        {
                            MPI_Send(&input[0], local_length, MPI_INT, rank, shift_start, comm);
                        }
                        else if (shift_end >= offsets_start[rank + 1])
                        {
                            MPI_Send(&input[0], offsets_end[rank] - shift_start + 1, MPI_INT, rank, shift_start, comm); //From calculated start to end of bucket -> Overlapping to next
                            MPI_Send(&input[offsets_end[rank] - shift_start + 1]
                                ,local_length - (offsets_end[rank] - shift_start + 1), MPI_INT, rank + 1, offsets_start[rank+1], comm); //MPI send end piece to rank where end is found; from start of bucket to calculated end
                            break;

                        }
                        else
                        {
                            std::cout << "Error 101 on rank "<<world_rank << std::endl;
                        }


                    }
                }
                if ((shift_start < 0) and (shift_end)>= 0)
                {
                        MPI_Send(&input[local_length - (shift_end + 1)], shift_end + 1, MPI_INT, 0, 0, comm);
                }


            }

            //Recieving
            if (((offsets_end[world_size - 1] - offsets_end[world_rank]) >= h)) {//No zeroes are needed

                if (local_length > h){//Local shifting from direct neghbour
                    MPI_Status status;
                    int position;
                    int number_amount;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);
                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start

                }
                else
                {
                    MPI_Status status;
                    int position;
                    int number_amount;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);
                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start

                    if (number_amount < local_length){ //If not enough recieve more
                        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        position = status.MPI_TAG;
                        MPI_Get_count(&status, MPI_INT, &number_amount);
                        MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start

                    }

                }
            }
            else if ((offsets_end[world_size -1] + 1 - offsets_start[world_rank] <= h))//Only zeroes are filled
            {
                std::fill(&input[local_length - h], &input[local_length], -1);
            }
            else { //Mixed case

                if (world_rank != world_size - 1){
                    MPI_Status status;
                    int position;
                    int number_amount;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);

                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start

                    //Fill the rest with zeroes
                    std::fill(&input[number_amount], &input[local_length], -1);
                }
                else
                {
                    std::fill(&input[local_length - h], &input[local_length], -1);
                }

            }

    }

// function to check if all singleton 
bool all_singleton (std::vector<tuple_ISA<int>> &input, MPI_Comm comm, const int world_rank,
             const int world_size, int local_length){

    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
    MPI_Datatype types_ISA[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths_ISA, displacements_ISA, types_ISA, &MPI_TUPLE_ISA);
    MPI_Type_commit(&MPI_TUPLE_ISA);

    if (world_rank != 0)
    {
        MPI_Send(&input[0], 1, MPI_TUPLE_ISA, world_rank - 1, 0, comm);
    }

    tuple_ISA<int> after_end[1];
    if (world_rank != world_size - 1)
    {
        MPI_Recv(after_end, 1, MPI_TUPLE_ISA, world_rank + 1, 0, comm, MPI_STATUS_IGNORE);
    }
    else
    {
        after_end[0].SA = -1;
    }


    bool singleton = true;

    for (int i = 0; i < local_length - 1; i++){
        if (input[i].B == input[i + 1].B)
            singleton = false;
    }
    if (input[local_length - 1].B == after_end[0].B)
        singleton = false;

    return singleton;
}
