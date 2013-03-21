#include <iostream>
#include <mpi.h>

using namespace std;

int main() {
  int sum = 0;
  int count = 0;
  int comm_sz;
  int my_rank;
  int block_size = 0;
  MPI_Comm comm;
  double start, finish;

  int *values = new int[100000];
  int *local_values = new int[1];

  for (int i = 0; i < 100000; ++i)
    values[i] = 0;

  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &my_rank);

  int *res = new int[comm_sz];
  for (int i = 0; i < comm_sz; ++i)
    res[i] = 0;

  if (my_rank != 0) {
    MPI_Barrier(comm);
    MPI_Recv(&block_size, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
    int *local_values = new int[block_size];
    MPI_Barrier(comm);
    MPI_Scatter(&(values[0]), block_size, MPI_INT, local_values, block_size, MPI_INT, 0, comm);
    int local_sum = 0;
    for (int i = 0; i < block_size; i++)
      local_sum += local_values[i];
    res[my_rank] = local_sum;
    MPI_Barrier(comm);
    MPI_Reduce(&(res[my_rank]), &sum, 1, MPI_INT, MPI_SUM, 0, comm);
  } else {
    int tmp;
    while (cin >> tmp) {
      values[count] = tmp;
      count++;
    }

    MPI_Barrier(comm);
    start = MPI_Wtime();
    block_size = count / comm_sz;

    for (int i = 0; i < comm_sz; ++i)
      MPI_Send(&block_size, 1, MPI_INT, i, 0, comm);
    MPI_Barrier(comm);
    local_values = new int[block_size];
    MPI_Scatter(&(values[0]), block_size, MPI_INT, &(local_values[0]), block_size, MPI_INT, 0, comm);
    
    for (int i = 0; i < block_size; ++i)
      res[my_rank] += local_values[i];

    MPI_Barrier(comm);
    int final_result = 0;
    MPI_Reduce(&(res[my_rank]), &final_result, 1, MPI_INT, MPI_SUM, 0, comm);

    sum += final_result;

    finish = MPI_Wtime();
    cout << "It took " << finish - start << " to get the sum " << sum << endl;
  }

  free(values);
  free(local_values);
  MPI_Finalize();

  return 0;
}
