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

tuple_vector get_kmers(const std::string& input, const int k){
    tuple_vector kmers;

    for (std::string::const_iterator i = input.begin(); i != input.end(); i += k) {
        kmers.push_back(
            std::make_tuple(
                std::string(i, i + k),
                int ((i - input.begin())/k)
            )
        );
    }
    return kmers;
}

void sort_tuple_vector(tuple_vector& input){
    std::sort(input.begin(), input.end(), [](std::tuple<std::string, int> a, std::tuple<std::string, int> b) {
        return std::lexicographical_compare(std::get<0>(a).begin(), std::get<0>(a).end(),
        std::get<0>(b).begin(), std::get<0>(b).end());
    });
}

void print_tuple_vector(const tuple_vector& input){

    for (tuple_vector::const_iterator i = input.begin(); i != input.end(); ++i) {
        std::cout << std::get<0>(*i) << "," << std::get<1>(*i) << "|";
    }
    std::cout << std::endl;
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