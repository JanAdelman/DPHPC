#include <unistd.h>
#include <sys/resource.h>

using namespace std;

// this file allows to track the memory usage of every node 
// example: showMemUsage("after open file ", world_rank); (shows the memory of each process after opening a file)

void showMemUsage(string s, int rank) {

  int tSize = 0, resident = 0, share = 0;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages                                                                                                                                                                                                                          
  double rss = resident * page_size_kb;
  double shared_mem = share * page_size_kb;
  std::cout << s << " on core " << rank << " RSS - " << rss/1000 << " Shared Memory - " << shared_mem/1000 << " Private Memory - "
  << (rss - shared_mem)/1000 << " MB\n";
}