//
// Chris Hansen
// Assignment 5
// Gaussian Elimination w/partial pivoting
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>


void create_matrix(double **matrix, double **matrix_copy);
void create_answer_vector(double *answer_vector);
void pivot_on_row(int row, double **matrix, double *answer_vector);
void convert_to_upper_triangle(int row, double **matrix, double *answer_vector);
void back_subsitution(double **matrix, double *answer_vector, double *answers);
int  find_largest_in_col(int col, double **matrix);

int SIZE_OF_MATRIX;
int NUMBER_OF_THREADS;
int BLOCK_SIZE;

int main(int argc, char *argv[]) {
  int i;
  int row;
  int col;

  if (argc != 3) {
    fprintf(stderr, "The arguments should be ./matrix_sum size_of_matrix number_of_threads\n");
    return 1;
  }

  SIZE_OF_MATRIX    = atoi(argv[1]);
  NUMBER_OF_THREADS = atoi(argv[2]);
  omp_set_num_threads(NUMBER_OF_THREADS);

  printf("Number of procs is %d\n",        omp_get_num_procs());
  printf("The number of threads is %d\n",  NUMBER_OF_THREADS);
  printf("Max number of threads is %d\n",  omp_get_max_threads());

  double **matrix = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double **matrix_copy = malloc(sizeof(double) *SIZE_OF_MATRIX);
  double *answer_vector = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answers = malloc(sizeof(double) * SIZE_OF_MATRIX);


  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    matrix[i] = malloc(sizeof(double) * SIZE_OF_MATRIX);
    matrix_copy[i] = malloc(sizeof(double) * SIZE_OF_MATRIX);
  }

  // Check to make sure we get an even block size
  if (SIZE_OF_MATRIX % NUMBER_OF_THREADS != 0) {
    fprintf(stderr, "The size of the matrix must be divisible by the number of threads\n");
    return 1;
  }

  BLOCK_SIZE = SIZE_OF_MATRIX / NUMBER_OF_THREADS;

  srand48(time(NULL)); // seed random number

  #pragma omp parallel num_threads(NUMBER_OF_THREADS) 
  {
    create_matrix(matrix, matrix_copy);
    create_answer_vector(answer_vector);
  }

  // Print matrix
  int j;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      printf("%g ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("Answers\n");

  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    printf("%g ", answer_vector[i]);
  }

  // Start Timing

  // Start elimination
  row = 0;
  col = 0;
  for (row = 0; row < SIZE_OF_MATRIX; ++row) {
    pivot_on_row(row, matrix, answer_vector);

    #pragma omp parallel num_threads(NUMBER_OF_THREADS)
      convert_to_upper_triangle(row, matrix, answer_vector);
  }
  // Print upper triangle matrix
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      printf("%g ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    printf("%g ", answer_vector[i]);
  }

  back_subsitution(matrix, answer_vector, answers);

  // Finish Timing

  double **l_norm = malloc(sizeof(double) * SIZE_OF_MATRIX);
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    l_norm[i] = malloc(sizeof(double) * SIZE_OF_MATRIX);
  }

  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      l_norm[i][j] = 0;
    }
  }

  int k;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      for (k = 0; k < SIZE_OF_MATRIX; ++k) {
        l_norm[i][j] += matrix[i][j] * answers[k];
      }
    }
  }

  double total = 0;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      total += l_norm[i][j] - matrix_copy[i][j];
    }
  }


  printf("L2 Norm is %g", total);

  //Ax - b = 0
  //x = the answer column at the end of upper triangle
  //b is the original matrix
  //A is the matrix after upper triangle

  printf("\n");
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    printf("%g\n", answers[i]);
  }
  printf("\n");

  // 
  // Cleanup
  //
  return 0;
}

void back_subsitution(double **matrix, double *answer_vector, double *answers) {
  double tmp_answer = 0;
  double tmp_subtraction = 0;
  int i;
  int j;

  for (i = SIZE_OF_MATRIX - 1; i >= 0 ; --i) {
    tmp_answer = answer_vector[i];
    // Subtraction
    for (j = SIZE_OF_MATRIX - 1; j > i; --j) {
      tmp_subtraction = matrix[i][j] * answers[j];
      tmp_answer = tmp_answer - tmp_subtraction;
    }

    // Division
    answers[i] = tmp_answer / matrix[i][i];
  }
}

void create_answer_vector(double *answer_vector) {
  int i;

  int my_thread_id = my_thread_id = omp_get_thread_num();

  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    double random = (drand48() - .5) * 2;
    //random = random * (double)10*10*10*10*10;
    answer_vector[i] = random;
  }
}

void create_matrix(double **matrix, double **matrix_copy) {
  int i;
  int j;
  int my_thread_id = omp_get_thread_num();
  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      double random = (drand48() - .5) * 2;
      //random = random * (double)10*10*10*10*10;
      matrix[i][j] = random;
      matrix_copy[i][j] = random;
      //matrix[i][j] = drand48();
    }
  }
}

void pivot_on_row(int row, double **matrix, double *answer_vector) { // CHECKING FOR LARGEST IS NOT WORKING RIGHT
  int i;
  int col = row; // This gives us 0,0 1,1 and so on
  //int *row_with_largest_in_col = malloc(sizeof(int) * NUMBER_OF_THREADS);
  int *largest_row = malloc(sizeof(int) * NUMBER_OF_THREADS);

  //
  // Find number of threads largest numbers with row
  //
  #pragma omp parallel num_threads(NUMBER_OF_THREADS)
    largest_row[omp_get_thread_num()] = find_largest_in_col(col, matrix);

  //
  // Find real largest row
  //

  int largest = 0;
  for (i = 0; i < NUMBER_OF_THREADS; ++i) {
    if (fabs(matrix[largest_row[i]][col]) > fabs(matrix[largest][col])) {
      largest = largest_row[i];
    }
  }

  // 
  // Pivot row with row_with_largest
  //
  if (row >= largest)
    return; // We're already the largest so no work needs to be done

  double *first_row = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *second_row = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answer_one = malloc(sizeof(double));
  double *answer_two = malloc(sizeof(double));
  answer_one[0] = answer_vector[row];
  answer_two[0] = answer_vector[largest];

  //
  // Copy rows
  //
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    first_row[i] = matrix[row][i];
    second_row[i] = matrix[largest][i];
  }

  //
  // Swap rows
  //
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    matrix[row][i] = second_row[i];
    matrix[largest][i] = first_row[i];
  }
  answer_vector[row] = answer_two[0];
  answer_vector[largest] = answer_one[0];
}

int find_largest_in_col(int col, double **matrix) {
  int i;
  int my_thread_id = omp_get_thread_num();
  double largest_num = 0;
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


void convert_to_upper_triangle(int row, double **matrix, double *answer_vector) {
  int i;
  int j;
  int col = row;
  int my_thread_id = omp_get_thread_num();
  double *row_copy = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answer_copy = malloc(sizeof(double));
  double *multipliers = malloc(sizeof(double) * SIZE_OF_MATRIX);

  // The row we're multiplying

  double denominator = matrix[row][col];
  double numerator;

  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    if (i <= row) {
      multipliers[i] = -1;
      continue;
    }

    numerator = matrix[i][col];
    multipliers[i] = numerator / denominator;

    for (j = 0; j < SIZE_OF_MATRIX; ++j)
      row_copy[j] = matrix[row][j];

    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      row_copy[j] = (row_copy[j] * numerator) / denominator;
    }

    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      matrix[i][j] = matrix[i][j] - row_copy[j];
    }
  }

  if (my_thread_id == 0) {
    for (i = row + 1; i < SIZE_OF_MATRIX; ++i) {
      if (multipliers[i] == -1) 
        continue;

      answer_copy[i] = answer_vector[row] * multipliers[i];
      answer_vector[i] = answer_vector[i] - answer_copy[i];
    }
  }

}

