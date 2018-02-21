#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <omp.h>
#include <math.h>
#include <assert.h>

/* Example filter sizes */
#define DATA_LEN  1024*512*512

/* Subtract the `struct timeval' values X and Y,
    storing the result in RESULT.
    Return 1 if the difference is negative, otherwise 0. */
int timeval_subtract (struct timeval * result, struct timeval * x, struct timeval * y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
    
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

double parallelFilterTest ( int data_len, unsigned int* input_array, unsigned int* output_array )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */
  #pragma omp parallel for 
  for (int x=0; x<data_len; x++) {
      /* it the data element matches the filter */ 
      if (input_array[x] > 1) {
      /* include it in the output */
      output_array[x] = input_array[x];
      }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Parallel took %lu seconds and %lu microseconds.\n", tresult.tv_sec, tresult.tv_usec );
  return dt;
}


int file_exist (char *filename)
{
  struct stat   buffer;   
  return (stat (filename, &buffer) == 0);
}

int main( int argc, char** argv )
{
  /* loop variables */
  int x,y;
  double dt;
  FILE *fp;

  /* Create matrixes */
  unsigned int * input_array;
  unsigned int * output_array;

  /* Initialize the data. Values don't matter much. */
  posix_memalign ( (void**)&input_array, 4096,  DATA_LEN * sizeof(unsigned int));
//  input_array = (unsigned int*) posix_memalign ( DATA_LEN * sizeof(unsigned int), 4096);
  for (x=0; x<DATA_LEN; x++)
  {
    input_array[x] = x % 2048;
  }

  /* start with an empty *all zeros* output array */
  posix_memalign ( (void**)&output_array, 4096, DATA_LEN * sizeof(unsigned int));
  memset ( output_array, 0, DATA_LEN );

  /* Create data directory */
  mkdir("data",0777);

  /* Loop Parallelism */
  printf("Begin estimating loop parallelism, data saved in data/test_parallel.txt\n");
  printf("\n");

  if (file_exist ("data/test_parallel.txt"))
    fp = fopen("data/test_parallel.txt","a+");
  else {
    fp = fopen("data/test_parallel.txt","w+");
    fprintf(fp,"num_threads,runtime\n");
  }

  for(int num_threads =1; num_threads<=16; num_threads*=2){
    printf("Num Threads: %d\n",num_threads);

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    fprintf(fp,"%d,",num_threads);
    
    dt = parallelFilterTest ( DATA_LEN, input_array, output_array);
    memset ( output_array, 0, DATA_LEN );
    fprintf(fp,"%f\n",dt);
  }
  fclose(fp);
}
