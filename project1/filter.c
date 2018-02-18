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

double DynamicFilter1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */
  #pragma omp parallel for schedule(dynamic) num_threads(16)
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

  printf ("Parallel(Dynamic) filter first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

double DynamicData1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
    /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the data */
  #pragma omp parallel for schedule(dynamic) num_threads(16)
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

  printf ("Parallel(Dynamic) data first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

double StaticFilter1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */
  #pragma omp parallel for schedule(static,128) num_threads(16)
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

  printf ("Parallel(Static) filter first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

double StaticData1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
    /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the data */
  #pragma omp parallel for schedule(static,128) num_threads(16)
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

  printf ("Parallel(Static) data first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

double LUFilter1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
  /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the filter */
  #pragma omp parallel for
  for (int y=0; y<filter_len; y+=4) { 
    for (int x=0; x<data_len; x+=4) {
      if (input_array[x] == filter_list[y])
        output_array[x] = input_array[x];
      if (input_array[x+1] == filter_list[y])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+2] == filter_list[y])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+3] == filter_list[y])
        output_array[x+3] = input_array[x+3];
    }
    for (int x=0; x<data_len; x+=4) {
      if (input_array[x] == filter_list[y+1])
        output_array[x] = input_array[x];
      if (input_array[x+1] == filter_list[y+1])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+2] == filter_list[y+1])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+3] == filter_list[y+1])
        output_array[x+3] = input_array[x+3];
    }
    for (int x=0; x<data_len; x+=4) {
      if (input_array[x] == filter_list[y+2])
        output_array[x] = input_array[x];
      if (input_array[x+1] == filter_list[y+2])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+2] == filter_list[y+2])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+3] == filter_list[y+2])
        output_array[x+3] = input_array[x+3];
    }
    for (int x=0; x<data_len; x+=4) {
      if (input_array[x] == filter_list[y+3])
        output_array[x] = input_array[x];
      if (input_array[x+1] == filter_list[y+3])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+2] == filter_list[y+3])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+3] == filter_list[y+3])
        output_array[x+3] = input_array[x+3];
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Parallel(Loop Unrolling) filter first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
  return dt;
}

double LUData1st ( int data_len, unsigned int* input_array, unsigned int* output_array, int filter_len, unsigned int* filter_list )
{
    /* Variables for timing */
  struct timeval ta, tb, tresult;
  double dt;

  /* get initial time */
  gettimeofday ( &ta, NULL );

  /* for all elements in the data */
  #pragma omp parallel for
  for (int x=0; x<data_len; x+=4) {
    for (int y=0; y<filter_len; y+=4) { 
      if (input_array[x] == filter_list[y])
        output_array[x] = input_array[x];
      if (input_array[x] == filter_list[y+1])
        output_array[x] = input_array[x];
      if (input_array[x] == filter_list[y+2])
        output_array[x] = input_array[x];
      if (input_array[x] == filter_list[y+3])
        output_array[x] = input_array[x];
    }
    for (int y=0; y<filter_len; y+=4) { 
      if (input_array[x+1] == filter_list[y])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+1] == filter_list[y+1])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+1] == filter_list[y+2])
        output_array[x+1] = input_array[x+1];
      if (input_array[x+1] == filter_list[y+3])
        output_array[x+1] = input_array[x+1];
    }
    for (int y=0; y<filter_len; y+=4) { 
      if (input_array[x+2] == filter_list[y])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+2] == filter_list[y+1])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+2] == filter_list[y+2])
        output_array[x+2] = input_array[x+2];
      if (input_array[x+2] == filter_list[y+3])
        output_array[x+2] = input_array[x+2];
    }
    for (int y=0; y<filter_len; y+=4) { 
      if (input_array[x+3] == filter_list[y])
        output_array[x+3] = input_array[x+3];
      if (input_array[x+3] == filter_list[y+1])
        output_array[x+3] = input_array[x+3];
      if (input_array[x+3] == filter_list[y+2])
        output_array[x+3] = input_array[x+3];
      if (input_array[x+3] == filter_list[y+3])
        output_array[x+3] = input_array[x+3];
    }
  }

  /* get initial time */
  gettimeofday ( &tb, NULL );

  timeval_subtract ( &tresult, &tb, &ta );
  dt = tresult.tv_sec + tresult.tv_usec / 1000000.0;  

  printf ("Parallel(Loop Unrolling) data first took %lu seconds and %lu microseconds.  Filter length = %d\n", tresult.tv_sec, tresult.tv_usec, filter_len );
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

  /* Create data directory */
  mkdir("data",0777);

  /* Execute at a variety of filter lengths */
  /* Loop Efficiency */
  printf("Begin estimating loop efficiency, data saved in data/loop_efficiency.txt\n");
  printf("\n");
  if (file_exist ("data/loop_efficiency.txt"))
    fp = fopen("data/loop_efficiency.txt","a+");
  else {
    fp = fopen("data/loop_efficiency.txt","w+");
    fprintf(fp,"filter_len,data_1st,filter_1st\n");
  }

  for ( int filter_len =1; filter_len<=FILTER_LEN; filter_len*=2) 
  {
    fprintf(fp, "%d,", filter_len);

    dt = serialDataFirst ( DATA_LEN, input_array, serial_array, filter_len, filter_list );
    fprintf(fp, "%f,", dt);
    memset ( output_array, 0, DATA_LEN );

    dt = serialFilterFirst ( DATA_LEN, input_array, output_array, filter_len, filter_list );
    checkData ( serial_array, output_array );
    fprintf(fp, "%f\n", dt);
    memset ( output_array, 0, DATA_LEN );
  }
  fclose(fp);

  /* Loop Parallelism */
  int fs = 512;
  printf("\n\nPreparing reference filtered data\n");
  memset ( serial_array, 0, DATA_LEN );
  serialDataFirst ( DATA_LEN, input_array, serial_array, fs, filter_list );

  printf("\n\n");
  printf("Begin estimating loop parallelism, data saved in data/loop_parallel.txt\n");
  printf("\n");

  if (file_exist ("data/loop_parallel.txt"))
    fp = fopen("data/loop_parallel.txt","a+");
  else {
    fp = fopen("data/loop_parallel.txt","w+");
    fprintf(fp,"num_threads,filter_1st,data_1st\n");
  }

  for(int num_threads =1; num_threads<=16; num_threads*=2){
    printf("Num Threads: %d\n",num_threads);

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    fprintf(fp,"%d,",num_threads);
    
    dt = parallelFilterFirst ( DATA_LEN, input_array, output_array, fs, filter_list );
    checkData ( serial_array, output_array );
    memset ( output_array, 0, DATA_LEN );
    fprintf(fp,"%f,",dt);

    dt = parallelDataFirst ( DATA_LEN, input_array, output_array, fs, filter_list );
    checkData ( serial_array, output_array );
    memset ( output_array, 0, DATA_LEN );
    fprintf(fp,"%f\n",dt);    
  }
  fclose(fp);

  /* Optimized Version */
  /* Loop unrolling */
  /* Here we run the function only on 16 threads with filter size 512 */
  printf("\n\n");
  printf("Begin estimating Loop Unrolling efficiency, data saved in data/loop_unrolling.txt\n");
  printf("\n");
  printf("Unrolling:4, Threads:16\n");

  if (file_exist ("data/loop_unrolling.txt"))
    fp = fopen("data/loop_unrolling.txt","a+");
  else {
    fp = fopen("data/loop_unrolling.txt","w+");
    fprintf(fp,"filter_1st,data_1st\n");
  }

  omp_set_dynamic(0);
  omp_set_num_threads(16);
  
  dt = LUFilter1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f,",dt);

  dt = LUData1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f\n",dt);    
  
  fclose(fp);

  /* Custom Scheduling */
  /* dynamic scheduling */
  printf("\n\n");
  printf("Begin estimating Custom Scheduling efficiency, data saved in data/dynamic_schedul.txt\n");
  printf("\n");
  printf("Threads:16\n");

  if (file_exist ("data/dynamic_schedul.txt"))
    fp = fopen("data/dynamic_schedul.txt","a+");
  else {
    fp = fopen("data/dynamic_schedul.txt","w+");
    fprintf(fp,"filter_1st,data_1st\n");
  }
  
  dt = DynamicFilter1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f,",dt);

  dt = DynamicData1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f\n",dt);    
  
  fclose(fp);

  /* Static scheduling with block size 128*/
  printf("\n\n");
  printf("Begin estimating Custom Scheduling efficiency, data saved in data/static_schedul.txt\n");
  printf("\n");
  printf("Static schedule, block size: 128\n");

  if (file_exist ("data/static_schedul.txt"))
    fp = fopen("data/static_schedul.txt","a+");
  else {
    fp = fopen("data/static_schedul.txt","w+");
    fprintf(fp,"filter_1st,data_1st\n");
  }
  
  dt = StaticFilter1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f,",dt);

  dt = StaticData1st ( DATA_LEN, input_array, output_array, fs, filter_list );
  checkData ( serial_array, output_array );
  memset ( output_array, 0, DATA_LEN );
  fprintf(fp,"%f\n",dt);    
  
  fclose(fp);
  
  printf("\n\n");

}

