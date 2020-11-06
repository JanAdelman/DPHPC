#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

#include "prefix_doubling.cpp"
#include "rebucket.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

int main (){
    std::string test = "bananaahahahhahahh";
    int k=2;
    tuple_vector test_vec =  get_kmers(test, k);
    print_tuple_vector(test_vec);
    sort_tuple_vector(test_vec);
    print_tuple_vector(test_vec);
    tuple_vector B = rebucket(test_vec);
    print_tuple_vector(B);   
    for(int h=k; h<=test.size();h*=2){
        tuple_vector B_new(B.size());
        for(int i=0;i<B.size();i++){
            B_new[std::get<1>(test_vec[i])]=B[i];
        }
        B=B_new;
        print_tuple_vector(B);
        triple_vector B2(B.size());
        for(int i=0;i<B.size();i++){
            if(i<B.size()-h){
                B2[i]= std::make_tuple(std::get<1>(B[i]),std::get<1>(B[i+h]),i);
            }
            else{
                 B2[i]=std::make_tuple(std::get<1>(B[i]),0,i);
            }
        }
    }


    return 0;
} 