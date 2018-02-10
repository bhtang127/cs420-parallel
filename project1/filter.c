/*******************************************************************************
*
*  Filter a large array based on the values in a second array.
*
********************************************************************************/

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
#define DATA_LEN  512*512*256
#define FILTER_LEN  1024


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

/* Function to apply the filter with the filter list in the outside loop */
double serialFilterFirst ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */ 
  for (int y=0; y<filter_len; y++) { 
    /* for all elements in the data */
    for (int x=0; x<data_len; x++) {
      /* it the data element matches the filter */ 
      if (input_array[x] == filter_list[y]) {
        /* include it in the output */
        output_array[x] = input_array[x];
      }
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;

  printf ("Serial filter first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}


/* Function to apply the filter with the filter list in the outside loop */
double serialDataFirst ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the data */
  for (int x=0; x<data_len; x++) {
    /* for all elements in the filter */ 
    for (int y=0; y<filter_len; y++) { 
      /* it the data element matches the filter */ 
      if (input_array[x] == filter_list[y]) {
        /* include it in the output */
        output_array[x] = input_array[x];
      }
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Serial data first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

/* Function to apply the filter with the filter list in the outside loop */
double parallelFilterFirst ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */
  #pragma omp parallel for 
  for (int y=0; y<filter_len; y++) { 
    /* for all elements in the data */
    for (int x=0; x<data_len; x++) {
      /* it the data element matches the filter */ 
      if (input_array[x] == filter_list[y]) {
        /* include it in the output */
        output_array[x] = input_array[x];
      }
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Parallel filter first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}


/* Function to apply the filter with the filter list in the outside loop */
double parallelDataFirst ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
    /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the data */
  #pragma omp parallel for
  for (int x=0; x<data_len; x++) {
    /* for all elements in the filter */ 
    for (int y=0; y<filter_len; y++) { 
      /* it the data element matches the filter */ 
      if (input_array[x] == filter_list[y]) {
        /* include it in the output */
        output_array[x] = input_array[x];
      }
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Parallel data first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}


void checkData ( unsigned int * serialarray, unsigned int * parallelarray )
{
  for (int i=0; i<DATA_LEN; i++)
  {
    if (serialarray[i] != parallelarray[i])
    {
      printf("Data check failed offset %d\n", i );
      return;
    }
  }
}


int main( int argc, char** argv )
{
  /* loop variables */
  int x,y;
  double dt;
  FILE *fp1, *fp2, *fp3, *fp4;

  /* Create matrixes */
  unsigned int * input_array;
  unsigned int * serial_array;
  unsigned int * output_array;
  unsigned int * filter_list;

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
  posix_memalign ( (void**)&serial_array, 4096, DATA_LEN * sizeof(unsigned int));
  memset ( serial_array, 0, DATA_LEN );

  /* Initialize the filter. Values don't matter much. */
  filter_list = (unsigned int*) malloc ( FILTER_LEN * sizeof(unsigned int));
  for (y=0; y<FILTER_LEN; y++)
  {
    filter_list[y] = y;
  }

  /* open file for saving data */
  fp1 = fopen("loop_efficiency.txt","w+");
  fp2 = fopen("loop_parallel.txt","w+");
  fprintf(fp1,"filter_len,\
               data_1,data_1_norm,\
               filter_1,filter_1_norm\n");
  fprintf(fp2,"num_threads,\
               pdata_1st,\
               pfilter_1st\n");

  // /* Execute at a variety of filter lengths */
  // /* Loop Efficiency */
  // for ( int filter_len =1; filter_len<=FILTER_LEN; filter_len*=2) 
  // {
  //   fprintf(fp1, "%d,", filter_len);

  //   dt = serialDataFirst ( DATA_LEN, input_array, serial_array, filter_len, filter_list );
  //   fprintf(fp1, "%f,", dt);
  //   fprintf(fp1, "%f,", dt / filter_len);
  //   memset ( output_array, 0, DATA_LEN );

  //   dt = serialFilterFirst ( DATA_LEN, input_array, output_array, filter_len, filter_list );
  //   checkData ( serial_array, output_array );
  //   fprintf(fp1, "%f,", dt);
  //   fprintf(fp1, "%f\n", dt / filter_len);
  //   memset ( output_array, 0, DATA_LEN );
  // }

  /* Loop Parallelism */
  int fs = 512;
  serialDataFirst ( DATA_LEN, input_array, serial_array, fs, filter_list );
  for(int num_threads =1; num_threads<=16; num_threads*=2){
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    fprintf(fp2,"%d,",num_threads);
    
    dt=0;
    for(int it=0; it<5; it++){
      dt += parallelFilterFirst ( DATA_LEN, input_array, output_array, fs, filter_list );
      checkData ( serial_array, output_array );
      memset ( output_array, 0, DATA_LEN );
    }
    fprintf(fp2,"%f,",dt/5);

    dt=0;
    for(int it=0; it<5; it++){
      dt += parallelDataFirst ( DATA_LEN, input_array, output_array, fs, filter_list );
      checkData ( serial_array, output_array );
      memset ( output_array, 0, DATA_LEN );
    }
    fprintf(fp2,"%f\n",dt/5);    
  }

  fclose(fp1);
  fclose(fp2);
}

