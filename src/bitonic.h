//http://monismith.info/cs599/examples/bitonicSort.c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/*
#define N 1024
#define X 200
*/

void bitonicSort(int, int *, int);
void bitonicMerge(int, int *, int);
void bitonicCompare(int, int *, int);
void swap(int *, int *);

/*
int main(int argc, char ** argv)
{
  int i;
  int * arr;

  srand(12345678);

  arr = (int *) malloc(sizeof(int)*N);
  for(i = 0; i < N; i++)
    arr[i] = rand()%X;

  bitonicSort(1, arr, N);

  for(i = 0; i < N; i++)
    printf("%d ", arr[i]);

  free((void *) arr);

  return 0;
}
*/

void bitonicSort(int up, int * sequence, int size)
{
  if(size == 1)
    return;
  else
  {
    bitonicSort(1, sequence,size/2);
    bitonicSort(0, sequence+size/2, size/2);
    bitonicMerge(up, sequence, size);
  }
}

void bitonicMerge(int up, int * sequence, int size)
{
  if(size == 1)
    return;
  else
  {
    bitonicCompare(up, sequence, size);
    bitonicMerge(up, sequence, size/2);
    bitonicMerge(up, sequence+size/2, size/2);
  }
}

void bitonicCompare(int up, int * sequence, int size)
{
  int distance = size/2;
  int * start, * end = sequence+distance;
  #pragma omp parallel for private(start) shared(sequence, end, up, distance)
  for(start = sequence; start < end; start++)
    if( (*start > *(start+distance)) == up)
      swap(start, start+distance);
}

void swap(int * x, int * y)
{
  int temp = *x;
  *x = *y;
  *y = temp;
}