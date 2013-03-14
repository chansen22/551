#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;

int main() {
  int sum = 0;
  int comm_sz;
  int my_rank;
  int block_size;
  int offset;
  MPI_Comm comm;
  double start, finish;

  vector<int> values;
  vector<int>::iterator iter;

  int *block = new int[1]; // Array of ints to be sent for work

  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &my_rank);

  int **blocks = new int*[comm_sz];  // Holds blocks, one for each process

  if (my_rank != 0) {
    MPI_Barrier(comm);
    MPI_Recv(&block_size, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
    int *block = new int[block_size+1];
    MPI_Recv(block, block_size, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);

    for (int block_number = 0; block_number < block_size; ++block_number) {
      sum += block[block_number];
    }


    MPI_Send(&sum, 1, MPI_INT, 0, 0, comm);

    MPI_Barrier(comm);
  } else {
    int tmp;
    while(cin >> tmp) {
      values.push_back(tmp);
    }

    MPI_Barrier(comm);
    start = MPI_Wtime();

    block_size = values.size() / comm_sz;

    for (int process_number = 0; process_number < comm_sz; ++process_number) {
      blocks[process_number] = new int[block_size];
      for (int block_number = 0; block_number < block_size; ++block_number) {
        offset = (process_number * block_size) + block_number;
        blocks[process_number][block_number] = values[offset];
      }
    }

    for (int process_number = 1; process_number < comm_sz; ++process_number) {
      MPI_Send(&block_size, 1, MPI_INT, process_number, 0, comm);
      MPI_Send(&(blocks[process_number][0]), block_size, MPI_INT, process_number, 0, comm);
    }

    for (int block_number = 0; block_number < block_size; ++block_number) {
      sum += blocks[0][block_number];
    }

    int proc_sum = 0;
    for (int process_number = 1; process_number < comm_sz; ++process_number) {
      MPI_Recv(&proc_sum, 1, MPI_INT, process_number, 0, comm, MPI_STATUS_IGNORE);
      sum += proc_sum;
    }

    finish = MPI_Wtime();
    cout << "It took " << finish - start << " to get the sum " << sum << endl;
    MPI_Barrier(comm);
  }

  delete(blocks);
  delete(block);
  MPI_Finalize();
  return 0;
}
