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


void pivot_on_row(int row, double **matrix, double *answer_vector);
void convert_to_upper_triangle(int row, double **matrix, double *answer_vector);
void back_subsitution(double **matrix, double *answer_vector, double *answers);

void fill_matrix(double **matrix, double **matrix_copy);
void fill_answer(double *answer_vector, double *answer_vector_copy);

int SIZE_OF_MATRIX;
int NUMBER_OF_THREADS;
int BLOCK_SIZE;

int main(int argc, char *argv[]) {
  int i;
  int row;
  int col;
  time_t start;
  time_t finish;

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

  double **matrix            = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double **matrix_copy       = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answer_vector      = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answer_vector_copy = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answers            = malloc(sizeof(double) * SIZE_OF_MATRIX);


  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    matrix[i]      = malloc(sizeof(double) * SIZE_OF_MATRIX);
    matrix_copy[i] = malloc(sizeof(double) * SIZE_OF_MATRIX);
  }

  srand48(time(NULL)); // seed random number

  fill_matrix(matrix, matrix_copy);
  fill_answer(answer_vector, answer_vector_copy);

  // Start Timing
  start = time(NULL);

  // Start elimination
  row = 0;
  col = 0;
  int j;
  for (row = 0; row < SIZE_OF_MATRIX; ++row) {
    pivot_on_row(row, matrix, answer_vector);

    #pragma omp parallel for
      for (i = 0; i < SIZE_OF_MATRIX; ++i) {
        convert_to_upper_triangle(row, matrix, answer_vector);
      }
  }
  back_subsitution(matrix, answer_vector, answers);

  // Finish Timing
  finish = time(NULL);

  double seconds = (double) difftime(finish, start);

  printf("Time Taken: %f\n", seconds);

  double l2 = 0;

  double total = 0;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      total = total + matrix_copy[i][j] * answers[j];
    }
    l2 = l2 + pow( (total - answer_vector_copy[i]), 2);
    total = 0;
  }

  l2 = sqrt(l2);

  printf("L2 norm is %g\n", l2);

  free(matrix);
  free(matrix_copy);
  free(answer_vector);
  free(answer_vector_copy);
  free(answers);

  return 0;
}

void fill_matrix(double **matrix, double ** matrix_copy) {
  int i, j;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      matrix[i][j] = (drand48() - .5) * 2 * 10*10*10*10*10;
      matrix_copy[i][j] = matrix[i][j];
    }
  }
}

void fill_answer(double *answer_vector, double *answer_vector_copy) {
  int i;
  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    answer_vector[i] = (drand48() - .5) * 2 * 10*10*10*10*10;
    answer_vector_copy[i] = answer_vector[i];
  }
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

void pivot_on_row(int row, double **matrix, double *answer_vector) { 
  int i;
  int col = row; // This gives us 0,0 1,1 and so on

  //
  // Find number of threads largest numbers with row
  //
  int largest = col;

  for (i = col; i < SIZE_OF_MATRIX; ++i) {
    if (fabs(matrix[i][col]) > fabs(matrix[largest][col])) {
      largest = i;
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

void convert_to_upper_triangle(int row, double **matrix, double *answer_vector) {
  int i;
  int j;
  int col = row;
  int my_thread_id    = omp_get_thread_num();
  double *row_copy    = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *answer_copy = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double *multipliers = malloc(sizeof(double) * SIZE_OF_MATRIX);

  // The row we're multiplying

  double denominator = matrix[row][col];
  double numerator;

  for (i = my_thread_id * BLOCK_SIZE; i < BLOCK_SIZE + BLOCK_SIZE * my_thread_id; ++i) {
    if (i <= row)
      continue;

    numerator = matrix[i][col];
    multipliers[i] = numerator / denominator;

    for (j = 0; j < SIZE_OF_MATRIX; ++j)
      row_copy[j] = matrix[row][j];

    for (j = 0; j < SIZE_OF_MATRIX; ++j)
      row_copy[j] = (row_copy[j] * numerator) / denominator;

    for (j = 0; j < SIZE_OF_MATRIX; ++j)
      matrix[i][j] = matrix[i][j] - row_copy[j];
  }

  for (i = row + 1; i < SIZE_OF_MATRIX; ++i) {
    answer_copy[i]   = answer_vector[row] * multipliers[i];
    answer_vector[i] = answer_vector[i] - answer_copy[i];
  }
}
