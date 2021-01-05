#include <iostream>
#include <map>
#include <tuple>
#include <bitset>

unsigned long int encode (char* input, int size){
    std::bitset<4> p;

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

    return p.to_ulong();
}

int main () {
    char* input = "CA";
    int size = 2;

    std::cout << encode(input, size) << std::endl;
}