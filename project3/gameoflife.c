// gameoflife.c
// Name: Bohao Tang
// JHED: btang11

//////////////////////////////////////////////////////////////////////////////
// The Rules
//
// For a space that is 'populated':
//     Each cell with one or no neighbors dies, as if by solitude. 
//     Each cell with four or more neighbors dies, as if by overpopulation. 
//     Each cell with two or three neighbors survives. 
// For a space that is 'empty' or 'unpopulated'
//     Each cell with three neighbors becomes populated.
/////////////////////////////////////////////////////////////////////////////// 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "mpi.h"

#define DEFAULT_ITERATIONS 64
#define GRID_WIDTH  256
#define DIM  16     // assume a square grid

int mod( int x, int m ) {
  // return x % m, but when x is negative
  // it will return the smallest non-negative l * m + x
  // therefore mod(-1, 5) will be 4
  if ( x >= 0 )
    return x % m;
  else {
    while ( x < 0 )
      x += m;
    return x;
  }
}

int adjacent_count ( int *prev_row, int *next_row, int *current_grid, int length, int dim, int index ) {
  // return the count of values in 8 adjacent point of given index 
  // in a grid of size length * dim
  // current_grid is a int sub-grid of size length * dim
  // prev_row and next_row is the two adjacent row of this sub-grid in the whole grid
  int row_index, col_index;
  int count;
  int *whole_grid;
  int width = (length + 2) * dim;

  whole_grid = (int *) malloc ( width * sizeof(int) );
  for(int i = 0; i < width; i++) {
    if ( i < dim )
      whole_grid[i] = prev_row[i];
    else if ( i < width - dim )
      whole_grid[i] = current_grid[i - dim];
    else
      whole_grid[i] = next_row[i - width + dim];
  }

  count = 0;
  row_index = index / dim;
  col_index = index - row_index * dim;
  row_index += 1;

  for ( int i = -1; i <= 1; i++)
    for ( int j = -1; j <= 1; j++)
      count += whole_grid[ (row_index + i) * dim + mod(col_index + j, dim) ];

  count -= whole_grid[row_index * dim + col_index];

  free(whole_grid);

  return count;
}

void game_of_life ( int *prev_row, int *next_row, int *current_grid, int length, int dim ) {
  // current_grid is a int sub-grid of size length * dim
  // prev_row and next_row is the two adjacent row of this sub-grid in the whole grid
  // this function update current_grid to be its next status in game of life

  int count;
  int width = length * dim;
  int *result;
  result = (int *) malloc( length * dim * sizeof(int) );
  
  for (int i = 0; i < width; i++) {
    count = adjacent_count ( prev_row, next_row, current_grid, length, dim, i );
    if ( current_grid[i] == 1 ) {
      if ( count < 2 || count > 3 )
        result[i] = 0;
      else
        result[i] = 1;
    }
    else {
      if ( count == 3 )
        result[i] = 1;
      else
        result[i] = 0;
    }
  }

  for ( int i=0; i < width; i++) {
    current_grid[i] = result[i];
  }
  free(result);

  return;
}

int main ( int argc, char** argv ) {

  int global_grid[ 256 ] = 
   {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

  // MPI Standard variable
  int num_procs;
  int ID, j;
  int iters = 0;
  int num_iterations;

  // Setup number of iterations
  if (argc == 1) {
    num_iterations = DEFAULT_ITERATIONS;
  }
  else if (argc == 2) {
    num_iterations = atoi(argv[1]);
  }
  else {
    printf("Usage: ./gameoflife <num_iterations>\n");
    exit(1);
  }

  // Messaging variables
  MPI_Status stat;
  // TODO add other variables as necessary
  int *local_grid;
  int *prev_row, *next_row;
  int *sent_first_row, *sent_last_row;
  int length, width;
  int prev, next;

  // MPI Setup
  if ( MPI_Init( &argc, &argv ) != MPI_SUCCESS )
  {
    printf ( "MPI_Init error\n" );
  }

  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs ); // Set the num_procs
  MPI_Comm_rank ( MPI_COMM_WORLD, &ID );

  assert ( DIM % num_procs == 0 );

  // TODO Setup your environment as necessary
  length = DIM / num_procs;
  width = length * DIM;
  next =  ( ID + 1 ) % num_procs;
  prev = ID == 0 ? num_procs -1 : ID-1;

  local_grid = (int *) malloc ( length * DIM * sizeof(int) );
  prev_row = (int *) malloc ( DIM * sizeof(int) );
  next_row = (int *) malloc ( DIM * sizeof(int) );
  sent_first_row = (int *) malloc ( DIM * sizeof(int) );
  sent_last_row = (int *) malloc ( DIM * sizeof(int) );



  for ( iters = 0; iters < num_iterations; iters++ ) {
    // TODO: Add Code here or a function call to you MPI code
    // Data patition and sending partition from 0 to other threads
    if ( num_procs > 1 ) {
      if ( ID == 0 ) {
        for ( int i = 0; i < width; i++) {
          local_grid[i] = global_grid[i];
        }
        for ( int id = 1; id < num_procs; id++ ) {
          MPI_Ssend ( global_grid + id * width, width, MPI_INT, id, 2, MPI_COMM_WORLD );
        }
      }
      else {
        MPI_Recv ( local_grid, width, MPI_INT, 0, 2, MPI_COMM_WORLD, &stat );
      }
    }
    else {
      // if num_procs = 1, no sending
      for ( int i = 0; i < width; i++) {
          local_grid[i] = global_grid[i];
        }
    }

    // Prepare to send the first and last row between processes
    for ( int i = 0; i < DIM; i++ ) {
      sent_first_row[i] = local_grid[i];
      sent_last_row[i] = local_grid[ i + width - DIM ];
    }

    // Send necessary message from even nodes to odd and then odd to even
    // Each local grid needs its prev and next row to update status
    if ( num_procs > 1 ) {
      if ( ID % 2 == 0 )
      {
        MPI_Ssend ( sent_first_row, DIM, MPI_INT, prev, 2, MPI_COMM_WORLD ); 
        MPI_Recv ( next_row, DIM, MPI_INT, next, 2, MPI_COMM_WORLD, &stat );
        MPI_Ssend ( sent_last_row, DIM, MPI_INT, next, 2, MPI_COMM_WORLD ); 
        MPI_Recv ( prev_row, DIM, MPI_INT, prev, 2, MPI_COMM_WORLD, &stat );
      }
      else 
      { 
        MPI_Recv ( next_row, DIM, MPI_INT, next, 2, MPI_COMM_WORLD, &stat );
        MPI_Ssend ( sent_first_row, DIM, MPI_INT, prev, 2, MPI_COMM_WORLD );
        MPI_Recv ( prev_row, DIM, MPI_INT, prev, 2, MPI_COMM_WORLD, &stat );
        MPI_Ssend ( sent_last_row, DIM, MPI_INT, next, 2, MPI_COMM_WORLD ); 
      }
    }
    else {
      // if num_procs = 1, no sending
      for ( int i = 0; i < DIM; i++ ) {
        prev_row[i] = sent_last_row[i];
        next_row[i] = sent_first_row[i];
      }
    }

    // Processing gameoflife seperatly
    game_of_life ( prev_row, next_row, local_grid, length, DIM );

    // Send local grid back to 0th process
    if ( num_procs > 1 ) {
      if ( ID != 0 ) {
        MPI_Ssend ( local_grid, width, MPI_INT, 0, 2, MPI_COMM_WORLD );
      }
      else {
        for ( int i = 0; i < width; i++ ) {
          global_grid[i] = local_grid[i];
        }
        for ( int id = 1; id < num_procs; id++ ) {
          MPI_Recv ( global_grid + id * width, width, MPI_INT, id, 2, MPI_COMM_WORLD, &stat );
        }
      }
    }
    else {
      // if num_procs = 1, no sending
      for ( int i = 0; i < width; i++ ) {
        global_grid[i] = local_grid[i];
      }
    }

    // Output the updated grid state
    if ( ID == 0 ) {
      printf ( "\nIteration %d: final grid:\n", iters );
      for ( j = 0; j < GRID_WIDTH; j++ ) {
        if ( j % DIM == 0 ) {
          printf( "\n" );
        }
        printf ( "%d  ", global_grid[j] );
      }
      printf( "\n" );
    }
  }

  // TODO: Clean up memory
  free(local_grid);
  free(prev_row);
  free(next_row);
  free(sent_first_row);
  free(sent_last_row);

  MPI_Finalize(); // finalize so I can exit
}
