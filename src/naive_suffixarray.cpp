//#define CATCH_CONFIG_ENABLE_BENCHMARKING
// Program that creates a suffix array from a string in a sequential manner
#include <iostream>
#include <cstring>
#include <algorithm>
//#include <benchmark/benchmark.h>

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
  // struct to store the suffixe and the indexes of length len
  struct suffix suffixes[len];

  // store all suffixes in an array such that we can retain its index when sorting
  // the suffixes

  for (int i = 0; i < len; ++i)
  {
    suffixes[i].index = i;
    suffixes[i].suffix_i = (str+i);
  }

  std::sort(suffixes, suffixes+len, cmp);

  // initialise new suffix array
  int *suffixArray = new int[len];
  for (int i = 0; i < len; i++)
      suffixArray[i] = suffixes[i].index;

  return  suffixArray;
}

// function to print an array
void printArray(int array[], int n)
{
  for (int i = 0; i < n; i++)
  {
    std::cout << array[i] << " ";
  }
  std::cout << std::endl;
}

void run() {
      char str[] = "bananaxanahanana$";
      int len = std::strlen(str);
      int *suffixArray = constructSuffixArray(str, len);
      
      std::cout << "suffix array for: " << str << std::endl;
      printArray(suffixArray, len);
      std::cout << "process finished" << std::endl;
}

int main(){
  run();
  return 0;
}

/*

static void BM_NaiveSuffixArray(benchmark::State& state) {
  // Perform setup here
  for (auto _ : state) {
  
      //This code is timed
      run();
  }
}

// Register the function as a benchmark
BENCHMARK(BM_NaiveSuffixArray);
// Run the benchmark
BENCHMARK_MAIN();

*/