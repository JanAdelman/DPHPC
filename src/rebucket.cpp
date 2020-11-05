#include <iostream>
#include <string>
#include <vector>
#include <tuple>

typedef std::vector<std::tuple<std::string,int>> tuple_vector;

tuple_vector rebucket(tuple_vector vec){
    tuple_vector B;
    std::string prev_string=std::get<0>(vec[0]);
    int prev_index=0;
    int count=0;
    for(tuple_vector::const_iterator i=vec.begin(); i!=vec.end();i++){
        if(std::get<0>(*i)==prev_string){
            B.push_back(std::make_tuple(std::get<0>(*i),prev_index));
        }
        else{
            B.push_back(std::make_tuple(std::get<0>(*i),count));
            prev_string=std::get<0>(*i);
            prev_index= count;
        }
        count++;
    }
    return B;
}