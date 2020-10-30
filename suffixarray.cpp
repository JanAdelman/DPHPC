// Program that creates a suffix array from a string in a sequential manner

#include <iostream>
#include <cstring>
#include <algorithm>

// The information of one suffix is stored in this structure
struct suffix
{
  int index;
  char *suffix_i;
};

// helper function to compare two suffixes, used by sort function
int cmp(struct suffix a, struct suffix b)
{
    return std::strcmp(a.suffix_i, b.suffix_i) < 0? 1 : 0;
}

// This function builds the suffix array. It returns the array with additional
// index
int *constructSuffixArray(char *str, int len)
{
  // struct to store the suffixex and the indexes of length len
  struct suffix suffixes[len];

  // store all suffixes in an array such that we can retain its index when sorting
  // the suffixes

  for (int i = 0; i < len; ++i)
  {
    suffixes[i].index = i;
    suffixes[i].suffix_i = (len+i);
  }

  std::sort(suffixes.begin(), suffixes.end(), cmp);


}


int main()
{
  char str[] = "hello";
  int len = std::strlen(str);
  std::cout << "process finished" << std::endl;
}
