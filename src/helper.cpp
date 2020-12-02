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

#define K 2

struct tuple_t
{
    int idx;
    char seq[K];
};

void print_char_array(const char *input, size_t size)
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

bool tuple_t_compare(const tuple_t &b, const tuple_t &a)
{
    return !std::lexicographical_compare(a.seq, a.seq + K, b.seq, b.seq + K);
}

void tuple_t_sort(tuple_t *input, size_t size)
{
    std::sort(input, input + size, tuple_t_compare);
}

void tuple_t_print(const tuple_t *input, size_t size)
{
    for (int i = 0; i < size; i++)
        print_char_array((input + i)->seq, K);
    std::cout << "---" << std::endl;
}

bool char_array_comp(char* a, char* b, int size){
    bool same = true;
    for (int i = 0; i < size; i++) {
        if (*(a + i) != *(b + i))
            same = false;
    }
    return same; 
}

// Adapted from http://selkie-macalester.org/csinparallel/modules/MPIProgramming/build/html/mergeSort/mergeSort.html
template <typename T>
tuple_t* typename_t_sort(int height, int id, T localArray[], int size, MPI_Comm comm)
{
    //SETUP TUPLE STRUCT
    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K};
    const MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    int parent, rightChild, myHeight;
    T *half1, *half2, *mergeResult;

    myHeight = 0;
    //PLEASE OVERLOAD

    half1 = localArray; // assign half1 to localArray

    while (myHeight < height)
    { // not yet at top
        parent = (id & (~(1 << myHeight)));

        if (parent == id)
        { // left child
            rightChild = (id | (1 << myHeight));

            int recieved_size;
            if (std::is_same<T, tuple_t>::value)
            {
                MPI_Status status;
                MPI_Probe(rightChild, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_TUPLE_STRUCT, &recieved_size);
                // allocate memory and receive array of right child
                half2 = (T *)malloc(recieved_size * sizeof(T));
                MPI_Recv(half2, recieved_size, MPI_TUPLE_STRUCT, rightChild, 0, MPI_COMM_WORLD, &status);
            }
            mergeResult = (T *)malloc((size+recieved_size)* sizeof(T));
            // merge half1 and half2 into mergeResult
            std::merge(half1, half1 + size, half2, half2 + recieved_size,mergeResult, tuple_t_compare);

            // reassign half1 to merge result
            half1 = mergeResult;

            free(half2);
            mergeResult = NULL;

            if (myHeight==1&&id==0)
            {

                return half1;
            }

            size=size+recieved_size;
            myHeight++;
        }
        else
        {   // right child
            // send local array to parent
            MPI_Send(half1, size, MPI_TUPLE_STRUCT, parent, 0, MPI_COMM_WORLD);
            if (myHeight != 0)
                free(half1);
            myHeight = height;
        }
    }
    return NULL;
}

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
            std::cout<<kmers[i].idx<<std::endl;
            local_SA[i]=kmers[i].idx;
        }
        else{
            index[i]=displ+i;
            std::cout<<kmers[i].idx<<std::endl;
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
