#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

/*
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch_amalgamated.hpp"
*/

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;
struct tuple_t
{
    int idx;
    char seq[];
};

/*
tuple_vector get_kmers(const std::string &input, const int k)
{
    tuple_vector kmers(input.size());

    for (std::string::const_iterator i = input.begin(); i != input.end(); i++)
    {
        kmers[i-input.begin()]=(
            std::make_tuple(
                std::string(i, i + k),
                int((i - input.begin()))));
    }
    return kmers;
}
*/

void get_kmers(const char* input, const int k, tuple_t* kmers)
{
    //tuple_t kmers[input.size()];

    for (int i = 0; i < sizeof(kmers); i++)
    {
        memcpy(kmers[i].seq, input + i, k);
        kmers[i].seq[strlen(kmers[i].seq)] = '\0'; /* Add terminator */
        
        kmers[i].idx = i;
    }
}

void print_char_array(const char* input){
    for (int i = 0; i < sizeof(input); ++i)
        std::cout << input[i]; 
    std::cout << std::endl;
}

void print_kmers(const tuple_t* input){
    for (int i = 0; i < sizeof(input); ++i)
        print_char_array(input[i].seq); 
    std::cout << "---" << std::endl; 
}

void sort_tuple_vector(tuple_vector &input)
{
    std::sort(input.begin(), input.end(), [&](auto a, auto b) {
        return std::lexicographical_compare(std::get<0>(a).begin(), std::get<0>(a).end(),
                                            std::get<0>(b).begin(), std::get<0>(b).end());
    });
}

void print_vector(const std::vector<int> &input){
    for(std::vector<int>::const_iterator i=input.begin();i!=input.end(); ++i){
        std::cout<<*i<<",";

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
    std::tuple<int, int, int> b){
    if(std::get<0>(a)<std::get<0>(b))
        return true;
    else if(std::get<0>(a)>std::get<0>(b))
        return false;
    else{
        if(std::get<1>(a)<std::get<1>(b))
            return true;
        else if(std::get<1>(a)>std::get<1>(b))
            return false;
        else{
            if(std::get<2>(a)<std::get<2>(b))
                return true;
            else// if(std::get<2>(a)>std::get<2>(b))
                return false;
        }
    }
}

void sort_triple_vector(triple_vector &input)
{   
    std::sort(input.begin(), input.end(), triple_tuple_comp);
    
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
            B[i-vec.begin()]=std::make_tuple(std::get<0>(*i), prev_index);
        }
        else
        {
            B[i-vec.begin()]=(std::make_tuple(std::get<0>(*i), count));
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
    int prev_index =0;
    int count = 0;
    for (triple_vector::const_iterator i = vec.begin(); i != vec.end(); i++)
    {
        if (std::get<0>(*i)!=prev_index_B or std::get<1>(*i)!=prev_index_B2)
        {
            B[i-vec.begin()]=count;
            prev_index_B=std::get<0>(*i);
            prev_index_B2=std::get<1>(*i);
            prev_index=count;
            SA[i-vec.begin()]=std::get<2>(*i);
        }
        else
        {
            B[i-vec.begin()]=prev_index;
            SA[i-vec.begin()]=std::get<2>(*i);
        }
        count++;
    }
}
bool check_singleton(std::vector<int> &input){
    for(std::vector<int>::const_iterator i=input.begin();i!=input.end()-1; ++i){
        if(*i==*(i+1)){
            return false;
        }
    }
    return true;
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