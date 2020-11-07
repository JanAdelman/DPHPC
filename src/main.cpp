#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

#include "helper.cpp"

typedef std::vector<std::tuple<std::string, int>> tuple_vector;
typedef std::vector<std::tuple<int, int, int>> triple_vector;

int main (){
    std::string test = "bananaxanahanana$";
    int k=5;
    tuple_vector test_vec =  get_kmers(test, k);
    //print_tuple_vector(test_vec);
    sort_tuple_vector(test_vec);
    //print_tuple_vector(test_vec);
    tuple_vector B_init = rebucket(test_vec);
    //print_tuple_vector(B);
    std::vector<int> SA(B_init.size());
    std::vector<int> B(B_init.size());   
    for(int h=k; h<=test.size();h*=2){
        std::vector<int> B_new(B.size());
        //reorder to string order
        for(int i=0;i<B.size();i++){
            if(h==k){
                //first iteration different because of datastructure
                B_new[std::get<1>(test_vec[i])]=std::get<1>(B_init[i]);       
            }
            else{
                //now we made SA a vector of ints
                B_new[SA[i]]=B[i];
            }
        }
        //Do we destroy B_new now? Consider std::move
        B=B_new;
        //print_tuple_vector(B);
        triple_vector BB2SA(B.size());
        for(int i=0;i<B.size();i++){
            if(i<B.size()-h){
                BB2SA[i]= std::make_tuple(B[i],B[i+h],i);
            }
            else{
                BB2SA[i]=std::make_tuple(B[i],0,i);
            }
        }
        print_triple_vector(BB2SA);
        sort_triple_vector(BB2SA);
        std::cout << "TRIPLE" << std::endl;
        print_triple_vector(BB2SA);
        rebucket_2h(BB2SA, SA, B);
        std::cout << "SA" << std::endl;
        print_vector(SA);
        std::cout << "B" << std::endl;
        print_vector(B);
        if(check_singleton(B)){
            print_vector(SA);
            return 0;
        }

    }

    return 0;
} 