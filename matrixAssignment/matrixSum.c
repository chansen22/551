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
void convert_to_upper_triangle(int row, int col, double **matrix, double *answer_vector);
void back_subsitution(double **matrix, double *answer_vector, double *answers);

void fill_matrix(double **matrix, double **matrix_copy);
void fill_answer(double *answer_vector, double *answer_vector_copy);

int SIZE_OF_MATRIX;
int NUMBER_OF_THREADS;
int BLOCK_SIZE;

int main(int argc, char *argv[]) {
  int i;
  int j;
  int row;

  double start, finish;

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
  //matrix[0][0] = 25;
  //matrix[0][1] = 5;
  //matrix[0][2] = 1;
  //matrix[1][0] = 64;
  //matrix[1][1] = 8;
  //matrix[1][2] = 1;
  //matrix[2][0] = 144;
  //matrix[2][1] = 12;
  //matrix[2][2] = 1;

  //answer_vector[0] = 106.8;
  //answer_vector[1] = 177.2;
  //answer_vector[2] = 279.2;

  //for (i = 0; i < 3; ++i)
  //  answer_vector_copy[i] = answer_vector[i];

  //for (i = 0; i < 3; ++i)
  //  for (j = 0; j < 3; ++j)
  //    matrix_copy[i][j] = matrix[i][j];

  // Start Timing
  start = omp_get_wtime();

  //printf("Starting\n");
  //for (i = 0; i < SIZE_OF_MATRIX; ++i) {
  //  for (j = 0; j < SIZE_OF_MATRIX + 1; ++j) {
  //    if (j == SIZE_OF_MATRIX)
  //      printf("|[%g]", answer_vector[i]);
  //    else
  //      printf("[%g] ", matrix[i][j]);
  //  }
  //  printf("\n");
  //}
  //printf("\n");

  // Start elimination
  row = 0;
  for (row = 0; row < SIZE_OF_MATRIX; ++row) {
    pivot_on_row(row, matrix, answer_vector);

    //printf("After Pivot\n");
    //for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    //  for (j = 0; j < SIZE_OF_MATRIX + 1; ++j) {
    //    if (j == SIZE_OF_MATRIX)
    //      printf("|[%g]", answer_vector[i]);
    //    else
    //      printf("[%g] ", matrix[i][j]);
    //  }
    //  printf("\n");
    //}
    //printf("\n");

    #pragma omp parallel for private(i)
      for (i = row + 1; i < SIZE_OF_MATRIX; ++i)
        convert_to_upper_triangle(i, row, matrix, answer_vector);

    //printf("After Annihilation\n");
    //for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    //  for (j = 0; j < SIZE_OF_MATRIX + 1; ++j) {
    //    if (j == SIZE_OF_MATRIX)
    //      printf("|[%g]", answer_vector[i]);
    //    else
    //      printf("[%g] ", matrix[i][j]);
    //  }
    //  printf("\n");
    //}
    //printf("\n");
    
  }
  back_subsitution(matrix, answer_vector, answers);

  //for (i = 0; i < SIZE_OF_MATRIX; ++i)
  //  printf("[%g] ", answers[i]);
  //printf("\n");

  // Finish Timing
  finish = omp_get_wtime();

  //for (i = 0; i < SIZE_OF_MATRIX; ++i) {
  //  for (j = 0; j < SIZE_OF_MATRIX; ++j) {
  //    printf("[%g] ", matrix[i][j]);
  //  }
  //  printf("\n");
  //}

  double seconds = finish - start;

  printf("Time Taken: %g\n", seconds);

  double total[SIZE_OF_MATRIX];

  double l2 = 0;

  for (i = 0; i < SIZE_OF_MATRIX; ++i) {
    total[i] = 0;
    for (j = 0; j < SIZE_OF_MATRIX; ++j) {
      total[i] += matrix_copy[i][j] * answers[j];
      //printf("Total at %d is %g\n", j, total[i]);
    }
    l2 += pow( (total[i] - answer_vector_copy[i]), 2);
    //printf("Answer vector copy is %g\n", answer_vector_copy[i]);
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
  int i;
  int j;

  for (i = SIZE_OF_MATRIX - 1; i >= 0 ; --i) {
    tmp_answer = answer_vector[i];
    // Subtraction
    for (j = SIZE_OF_MATRIX - 1; j > i; --j)
      tmp_answer = tmp_answer - (matrix[i][j] * answers[j]);

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

  double *tmp;
  double tmpAnswer;

  tmp = matrix[row];
  matrix[row] = matrix[largest];
  matrix[largest] = tmp;
  tmpAnswer = answer_vector[largest];
  answer_vector[largest] = answer_vector[row];
  answer_vector[row] = tmpAnswer;

}

void convert_to_upper_triangle(int row, int col, double **matrix, double *answer_vector) {
  int j;
  double *row_copy = malloc(sizeof(double) * SIZE_OF_MATRIX);
  double answer_copy;

  double denominator = matrix[col][col];
  double numerator = matrix[row][col];

  for (j = 0; j < SIZE_OF_MATRIX; ++j) {
    row_copy[j] = matrix[col][j] * numerator / denominator;
    matrix[row][j] = matrix[row][j] - row_copy[j];
  }

  answer_copy = answer_vector[col] * numerator / denominator;
  answer_vector[row] = answer_vector[row] - answer_copy;

  free(row_copy);
}
