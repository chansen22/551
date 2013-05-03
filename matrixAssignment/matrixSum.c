//
// Chris Hansen
// Assignment 5
// Gaussian Elimination w/partial pivoting
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


void create_matrix(double **matrix);
void pivot_on_row(int row, double **matrix);
void convert_to_upper_triangle(int row, double **matrix);
int  find_largest_in_col(int col, double **matrix);

int SIZE_OF_MATRIX;
int NUMBER_OF_THREADS;
int BLOCK_SIZE;

int main(int argc, char *argv[]) {
  int i;
  int row;
  int col;

  if (argc != 3) {
    fprintf(stderr, "The arguments should be ./matrix_sum size_of_matrix NUMBER_OF_THREADS\n");
    return 1;
  }

  SIZE_OF_MATRIX    = atoi(argv[1]);
  NUMBER_OF_THREADS = atoi(argv[2]);

  printf("The size of the matrix is %d\n", SIZE_OF_MATRIX);
  printf("The number of threads is %d\n",  NUMBER_OF_THREADS);
  printf("Max number of threads is %d\n",  omp_get_max_threads());

  double **matrix = malloc(sizeof(double) * SIZE_OF_MATRIX);
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    matrix[i] = malloc(sizeof(double) * SIZE_OF_MATRIX);
  }

  double *variables = malloc(sizeof(double) * SIZE_OF_MATRIX);

  // Check to make sure we get an even block size
  if (SIZE_OF_MATRIX % NUMBER_OF_THREADS != 0) {
    fprintf(stderr, "The size of the matrix must be divisible by the number of threads\n");
    return 1;
  }

  BLOCK_SIZE = SIZE_OF_MATRIX / NUMBER_OF_THREADS;

  srand48(42); // seed random number

  #pragma omp parallel num_threads(NUMBER_OF_THREADS)
    create_matrix(matrix);

  // Variable matrix
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    variables[i] = ((drand48() - .5) * 2) * 10*10*10*10*10*10*10*10*10*10;
  }

  // Print matrix
  int j;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      printf("[%g]", matrix[i][j]);
    }
    printf("\n");
  }

  // Start elimination
  row = 0;
  col = 0;
  for (row = 0; row < SIZE_OF_MATRIX; ++row) {
    pivot_on_row(row, matrix);

    #pragma omp parallel num_threads(NUMBER_OF_THREADS)
      convert_to_upper_triangle(row, matrix);
      // split into blocks (if) and multiply the original row, i, by work row
      // subtract and then save in original matrix
  }

  printf("After swap\n");
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      printf("[%g]", matrix[i][j]);
    }
    printf("\n");
  }

  // 
  // Cleanup
  //
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    free(matrix[i]);
  }
  free(matrix);
  free(variables);
  return 0;
}

void create_matrix(double **matrix) {
  int i;
  int j;
  int my_thread_id = omp_get_thread_num();
  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      double random = (drand48() - .5) * 2;
      random = random * (double)10*10*10*10*10*10*10*10*10*10;
      //matrix[i][j] = random;
      matrix[i][j] = drand48();
    }
  }
}

void pivot_on_row(int row, double **matrix) {
  int i;
  int col = row; // This gives us 0,0 1,1 and so on
  int *row_with_largest_in_col = malloc(sizeof(int) * NUMBER_OF_THREADS);

  //
  // Find number of threads largest numbers with row
  //
  #pragma omp parallel num_threads(NUMBER_OF_THREADS)
    row_with_largest_in_col[omp_get_thread_num()] = find_largest_in_col(col, matrix);

  //
  // Find real largest row
  //
  int row_with_largest = 0;
  for (i = 0; i < NUMBER_OF_THREADS; ++i) {
    if (row_with_largest_in_col[i] >= row_with_largest_in_col[row_with_largest])
      row_with_largest = row_with_largest_in_col[i];
  }

  printf("Swapping row %d with row %d\n", row, row_with_largest);

  // 
  // Pivot row with row_with_largest
  //
  if (row == row_with_largest)
    return; // We're already the largest so no work needs to be done

  double *first_row = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *second_row = malloc(sizeof(double) * SIZE_OF_MATRIX);

  //
  // Copy rows
  //
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    first_row[i] = matrix[row][i];
    second_row[i] = matrix[row_with_largest][i];
  }

  //
  // Swap rows
  //
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    matrix[row][i] = second_row[i];
    matrix[row_with_largest][i] = first_row[i];
  }
  
  // Cleanup
  free(first_row);
  free(second_row);
  free(row_with_largest_in_col);
}

int find_largest_in_col(int col, double **matrix) {
  int i;
  int my_thread_id = omp_get_thread_num();
  double largest_num = fabs(matrix[col][col]);
  int return_row = 0;

  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    if (!(i < col)) { // This is a hack to not allow the swap to swap above what we're working on it wastes threads
      if (fabs(matrix[i][col]) >= largest_num) {
        largest_num = fabs(matrix[i][col]);
        return_row = i;
      }
    }
  }
  return return_row;
}


void convert_to_upper_triangle(int row, double **matrix) {
  printf("Converting to upper triangle, row is %d\n", row);
  int i;
  int col = row;
  int my_thread_id = omp_get_thread_num();
  double *row_copy = malloc(sizeof(double) * SIZE_OF_MATRIX);

  for (i = 0; i < SIZE_OF_MATRIX; ++i)
    row_copy[i] = matrix[row][i];

  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    
  }

  free(row_copy);
}



