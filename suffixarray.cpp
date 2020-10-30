// Program that creates a suffix array from a string in a sequential manner

#include <iostream>
#include <cstring>
#include <algorithm>

// The information of one suffix is stored in this structure
struct suffix
{
  int index;
  char *suffix;
};

// This function builds the suffix array. It returns the array with additional
// index
int constructSuffixArray(char *str, int len)
{
  
}


int main()
{
  char str[] = "hello";
  int len = std::strlen(str);
  std::cout << "process finished" << std::endl;
}
