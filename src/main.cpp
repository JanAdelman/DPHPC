#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

#include "prefix_doubling.cpp"
#include "rebucket.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;

int main (){
    std::string test = "bananaahahahhahahh";

    tuple_vector test_vec =  get_kmers(test, 2);
    print_tuple_vector(test_vec);
    sort_tuple_vector(test_vec);
    print_tuple_vector(test_vec);
    tuple_vector B = rebucket(test_vec);
    print_tuple_vector(B);
} 