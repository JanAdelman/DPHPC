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

tuple_vector get_kmers(const std::string &input, const int k)
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

void sort_tuple_vector(tuple_vector &input)
{
    std::sort(input.begin(), input.end(), [&](auto a, auto b) {
        return std::lexicographical_compare(std::get<0>(a).begin(), std::get<0>(a).end(),
                                            std::get<0>(b).begin(), std::get<0>(b).end());
    });
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
            else if(std::get<2>(a)>std::get<2>(b))
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
    tuple_vector B;
    std::string prev_string = std::get<0>(vec[0]);
    //Maybe not keep track of prev_index, but access i-1.
    int prev_index = 0;
    int count = 0;
    
    for (tuple_vector::const_iterator i = vec.begin(); i != vec.end(); i++)
    {
        if (std::get<0>(*i) == prev_string)
        {
            //Not use push_back, but allocate size
            B.push_back(std::make_tuple(std::get<0>(*i), prev_index));
        }
        else
        {
            B.push_back(std::make_tuple(std::get<0>(*i), count));
            prev_string = std::get<0>(*i);
            prev_index = count;
        }
        count++;
    }
    return B;
}

triple_vector rebucket_2h(triple_vector vec)
{
    //problem:Need B to be tuple vector in order to return at the
    //beginning of for loop, but we only saved the indices here in the 
    //triple vectors
    triple_vector B;
    std::string prev_string = std::get<0>(vec[0]);
    int prev_index = 0;
    int count = 0;
    for (tuple_vector::const_iterator i = vec.begin(); i != vec.end(); i++)
    {
        if (std::get<0>(*i) == prev_string)
        {
            B.push_back(std::make_tuple(std::get<0>(*i), prev_index));
        }
        else
        {
            B.push_back(std::make_tuple(std::get<0>(*i), count));
            prev_string = std::get<0>(*i);
            prev_index = count;
        }
        count++;
    }
    return B;
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