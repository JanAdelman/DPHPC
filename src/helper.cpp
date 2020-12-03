#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>
#include <type_traits>

/*
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch_amalgamated.hpp"
*/

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

#define K 1

struct tuple_t
{
    int idx;
    char seq[K];
};

struct triple_t
{
    int idx;
    int b;
    int b2;
};

void rebucketing(int *index,tuple_t *kmers, size_t size, int displ, int *local_SA){
    index[0]=displ;
    local_SA[0]=kmers[0].idx;
    for(int i=1;i<size;i++){
        if(!std::lexicographical_compare(kmers[i-1].seq, kmers[i-1].seq + K, kmers[i].seq, kmers[i].seq + K)){
            index[i]=index[i-1];
            //std::cout<<kmers[i].idx<<std::endl;
            local_SA[i]=kmers[i].idx;
        }
        else{
            index[i]=displ+i;
            //std::cout<<kmers[i].idx<<std::endl;
            local_SA[i]=kmers[i].idx;
        }
    }
}

void print_char_array(const char *input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout << input[i];
    std::cout << std::endl;
}

void print_int_array(const int *input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout << input[i];
    std::cout << std::endl;
}

void get_kmers(const char *input, const int k, tuple_t *kmers, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        memcpy(kmers[i].seq, input + i, k);
        kmers[i].seq[strlen(kmers[i].seq)] = '\0'; /* Add terminator */
        kmers[i].idx = i;
    }
}

void get_kmers_adapt(const char *input, const int k, tuple_t *kmers, size_t size, int displacement)
{
    for (int i = 0; i < size; i++)
    {
        memcpy(kmers[i].seq, input + i, k);
        kmers[i].seq[strlen(kmers[i].seq)] = '\0'; /* Add terminator */
        kmers[i].idx =displacement+ i;
    }
}

bool char_array_comp(const char* a, const char* b, int size){
    bool same = true;
    for (int i = 0; i < size; i++) {
        if (*(a + i) != *(b + i))
            same = false;
    }
    return same; 
}

bool tuple_t_compare(const tuple_t &a, const tuple_t &b)
{
    
    if(std::lexicographical_compare(a.seq, a.seq + K, b.seq, b.seq + K)){
        return true;
    }
    else{
        if(char_array_comp(a.seq,b.seq,K)){
            return a.idx<b.idx;
        }
        else{
            return false;
        }
    }
    
    //return !std::lexicographical_compare(b.seq, b.seq + K, a.seq, a.seq + K);
}
bool triple_t_compare(const triple_t &b, const triple_t &a)
{
    if (b.b < a.b)
        return true;
    else if (a.b > b.b)
        return false;
    else
    {
        if (a.b2 < b.b2)
            return true;
        else if (a.b2 > b.b2)
            return false;
        else
        {
            if (a.idx < b.idx)
                return true;
            else // if(std::get<2>(a)>std::get<2>(b))
                return false;
        }
    }
}

void t_sort(tuple_t *input, size_t size)
{
    std::sort(input, input + size, tuple_t_compare);
}
void t_sort(triple_t *input, size_t size)
{
    std::sort(input, input + size, triple_t_compare);
}

void t_print(const tuple_t *input, size_t size)
{
    for (int i = 0; i < size; i++){
        print_char_array((input + i)->seq, K);
        std::cout<<input[i].idx<<std::endl;
    }
    std::cout << "---" << std::endl;
}
void t_print(const triple_t *input, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        std::cout << "B: " << (input + i)->b << std::endl;
        std::cout << "B2: " << (input + i)->b2 << std::endl;
        std::cout << "SA: " << (input + i)->idx << std::endl;
    }
    std::cout << "---" << std::endl;
}



// Adapted from http://selkie-macalester.org/csinparallel/modules/MPIProgramming/build/html/mergeSort/mergeSort.html
tuple_t *typename_t_sort(int height, int id, tuple_t localArray[], int size, MPI_Comm comm)
{
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K};
    const MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    int parent, rightChild, local_height;
    tuple_t *half1, *half2, *mergeResult;

    local_height = 0;

    half1 = localArray; // assign half1 to localArray

    while (local_height < height)
    { // not yet at top
        parent = (id & (~(1 << local_height)));

        if (parent == id)
        { // left child
            rightChild = (id | (1 << local_height));

            int recieved_size;
            MPI_Status status;
            MPI_Probe(rightChild, 0, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_TUPLE_STRUCT, &recieved_size);
            // allocate memory and receive array of right child
            half2 = (tuple_t *)malloc(recieved_size * sizeof(tuple_t));
            MPI_Recv(half2, recieved_size, MPI_TUPLE_STRUCT, rightChild, 0, MPI_COMM_WORLD, &status);

            mergeResult = (tuple_t *)malloc((size + recieved_size) * sizeof(tuple_t));
            // merge half1 and half2 into mergeResult

            std::merge(half1, half1 + size, half2, half2 + recieved_size, mergeResult, tuple_t_compare);

            // reassign half1 to merge result
            half1 = mergeResult;

            free(half2);
            mergeResult = NULL;

            if (local_height == 1 && id == 0)
            {
                return half1;
            }

            size = size + recieved_size;
            local_height++;
        }
        else
        { // right child
            // send local array to parent
            MPI_Send(half1, size, MPI_TUPLE_STRUCT, parent, 0, MPI_COMM_WORLD);

            if (local_height != 0)
                free(half1);
            local_height = height;
        }
    }
    return NULL;
}
triple_t *typename_t_sort(int height, int id, triple_t localArray[], int size, MPI_Comm comm)
{
    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    const MPI_Aint displacements_triple[3] = {0, sizeof(int), sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);

    int parent, rightChild, local_height;
    triple_t *half1, *half2, *mergeResult;

    local_height = 0;

    half1 = localArray; // assign half1 to localArray

    while (local_height < height)
    { // not yet at top
        parent = (id & (~(1 << local_height)));

        if (parent == id)
        { // left child
            rightChild = (id | (1 << local_height));

            int recieved_size;
            MPI_Status status;
            MPI_Probe(rightChild, 0, comm, &status);

            MPI_Get_count(&status, MPI_TRIPLE_STRUCT, &recieved_size);
            // allocate memory and receive array of right child
            half2 = (triple_t *)malloc(recieved_size * sizeof(triple_t));
            MPI_Recv(half2, recieved_size, MPI_TRIPLE_STRUCT, rightChild, 0, comm, &status);

            mergeResult = (triple_t *)malloc((size + recieved_size) * sizeof(triple_t));
            // merge half1 and half2 into mergeResult

            std::merge(half1, half1 + size, half2, half2 + recieved_size, mergeResult, triple_t_compare);

            // reassign half1 to merge result
            half1 = mergeResult;

            free(half2);
            mergeResult = NULL;

            if (local_height == 1 && id == 0)
            {
                return half1;
            }

            size = size + recieved_size;
            local_height++;
        }
        else
        { // right child
            // send local array to parent
            MPI_Send(half1, size, MPI_TRIPLE_STRUCT, parent, 0, comm);

            if (local_height != 0)
                free(half1);
            local_height = height;
        }
    }
    return NULL;
}

void shift_h(int *input, const int h, MPI_Comm comm, const int world_rank,
             const int world_size, int *offsets, int local_length)
{

    int local_end_len = offsets[world_rank];
    int deleted_size = 0;

    for (int i = 0; i < world_size; i++){
        if (offsets[i] < h) //This offset will still be deleted
            deleted_size = offsets[i] + 1;
    }
   if (h > local_end_len) // This section is entirely shifted away
   { 
        std::fill(input,
            input + local_length, 0);
   }
    else {
        int shift_by = h - deleted_size;

        if ((world_rank != 0) and (offsets[world_rank - 1]>=h))
        {
            MPI_Send(input, shift_by, MPI_INT, world_rank - 1, 0, comm);
        }

        std::copy(input + shift_by, input + local_length, input);

        if (world_rank != world_size - 1)
        {
            MPI_Recv(input + local_length - shift_by, shift_by,
                     MPI_INT,
                     world_rank + 1, 0,
                     comm, MPI_STATUS_IGNORE);
            //Write into input at end
        }
        else
        {
            std::fill(input + local_length - shift_by,
                      input + local_length, 0);
        }
    }
    print_int_array(input, local_length);
}

bool all_singleton (int *input, MPI_Comm comm, const int world_rank,
             const int world_size, int local_length){

    if (world_rank != 0)
    {
        MPI_Send(input, 1, MPI_INT, world_rank - 1, 0, comm);
    }

    int after_end[1] = {-1}; 
    if (world_rank != world_size - 1)
    {
        MPI_Recv(after_end, 1, MPI_INT, world_rank + 1, 0, comm, MPI_STATUS_IGNORE);
    }

    bool singleton = true;

    for (int i = 0; i < local_length - 1; i++){
        if (input[i] == input[i + 1])
            singleton = false;
    }
    if (input[local_length - 1] == after_end[0])
        singleton = false;

    std::cout << singleton <<std::endl;
    return singleton; 
}

/*
void sort_tuple_vector(tuple_vector &input)
{
    std::sort(input.begin(), input.end(), [&](auto a, auto b) {
        return std::lexicographical_compare(std::get<0>(a).begin(), std::get<0>(a).end(),
                                            std::get<0>(b).begin(), std::get<0>(b).end());
    });
}

void print_vector(const std::vector<int> &input)
{
    for (std::vector<int>::const_iterator i = input.begin(); i != input.end(); ++i)
    {
        std::cout << *i << ",";
    }
    std::cout << std::endl;
}

void print_tuple_vector(const tuple_vector &input)
{

    for (tuple_vector::const_iterator i = input.begin(); i != input.end(); ++i)
    {
        std::cout << std::get<0>(*i) << "," << std::get<1>(*i) << "|";
    }
    std::cout << std::endl;
}

void print_triple_vector(const triple_vector &input)
{

    for (triple_vector::const_iterator i = input.begin(); i != input.end(); ++i)
    {
        std::cout << std::get<0>(*i) << "," << std::get<1>(*i) << "," << std::get<2>(*i) << "|";
    }
    std::cout << std::endl;
}

bool triple_tuple_comp(std::tuple<int, int, int> a,
                       std::tuple<int, int, int> b)
{
    if (std::get<0>(a) < std::get<0>(b))
        return true;
    else if (std::get<0>(a) > std::get<0>(b))
        return false;
    else
    {
        if (std::get<1>(a) < std::get<1>(b))
            return true;
        else if (std::get<1>(a) > std::get<1>(b))
            return false;
        else
        {
            if (std::get<2>(a) < std::get<2>(b))
                return true;
            else // if(std::get<2>(a)>std::get<2>(b))
                return false;
        }
    }
}

void sort_triple_vector(triple_vector &input)
{
    std::sort(input.begin(), input.end(), triple_tuple_comp);
}

void rebucketing(int *index,tuple_t *kmers, size_t size, int displ, int *local_SA){
    index[0]=displ;
    local_SA[0]=kmers[0].idx;
    for(int i=1;i<size;i++){
        if(!std::lexicographical_compare(kmers[i-1].seq, kmers[i-1].seq + K, kmers[i].seq, kmers[i].seq + K)){
            index[i]=index[i-1];
            //std::cout<<kmers[i].idx<<std::endl;
            local_SA[i]=kmers[i].idx;
        }
        else{
            index[i]=displ+i;
            //std::cout<<kmers[i].idx<<std::endl;
            local_SA[i]=kmers[i].idx;
        }
    }
}

tuple_vector rebucket(tuple_vector vec)
{
    std::string prev_string = std::get<0>(vec[0]);
    //Maybe not keep track of prev_index, but access i-1.
    int prev_index = 0;
    int count = 0;
    tuple_vector B(vec.size());

    for (tuple_vector::const_iterator i = vec.begin(); i != vec.end(); i++)
    {
        if (std::get<0>(*i) == prev_string)
        {
            //Not use push_back, but allocate size
            B[i - vec.begin()] = std::make_tuple(std::get<0>(*i), prev_index);
        }
        else
        {
            B[i - vec.begin()] = (std::make_tuple(std::get<0>(*i), count));
            prev_string = std::get<0>(*i);
            prev_index = count;
        }
        count++;
    }
    return B;
}

void rebucket_2h(triple_vector vec, std::vector<int> &SA, std::vector<int> &B)
{
    //problem:Need B to be tuple vector in order to return at the
    //beginning of for loop, but we only saved the indices here in the
    //triple vectors
    int prev_index_B = 0;
    int prev_index_B2 = 0;
    int prev_index = 0;
    int count = 0;
    for (triple_vector::const_iterator i = vec.begin(); i != vec.end(); i++)
    {
        if (std::get<0>(*i) != prev_index_B or std::get<1>(*i) != prev_index_B2)
        {
            B[i - vec.begin()] = count;
            prev_index_B = std::get<0>(*i);
            prev_index_B2 = std::get<1>(*i);
            prev_index = count;
            SA[i - vec.begin()] = std::get<2>(*i);
        }
        else
        {
            B[i - vec.begin()] = prev_index;
            SA[i - vec.begin()] = std::get<2>(*i);
        }
        count++;
    }
}
bool check_singleton(std::vector<int> &input)
{
    for (std::vector<int>::const_iterator i = input.begin(); i != input.end() - 1; ++i)
    {
        if (*i == *(i + 1))
        {
            return false;
        }
    }
    return true;
}

tuple_vector get_kmers_not_mpi(const std::string &input, const int k)
 {
     tuple_vector kmers;

     for (std::string::const_iterator i = input.begin(); i != input.end(); i++)
     {
         kmers.push_back(
             std::make_tuple(
                 std::string(i, i + k),
                 int((i - input.begin()))));
     }
     return kmers;
 }


*/

/*
int main (){
    std::string test = "bananaahahahhahahh";
    print_tuple_vector(get_kmers(test, 2));
}

std::tuple<std::string, int> data[] = {std::make_tuple("ba",0),std::make_tuple("an",1),std::make_tuple("na",1)};
tuple_vector test (data, data + sizeof(data) / sizeof(int) );

TEST_CASE( "Kmers are computed!", "[kmer]" ) {
    REQUIRE( get_kmers("bana", 2) == test );
}

*/
