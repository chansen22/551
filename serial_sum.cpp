#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

int main() {
  int sum = 0;
  int tmp;
  double start, finish;
  vector<int> values;
  vector<int>::iterator it;
  while(cin >> tmp) {
    values.push_back(tmp);
  }

  // Start time
  start = MPI_Wtime();

  for(it = values.begin(); it != values.end(); ++it){
    sum = sum + *it;
  }

  finish = MPI_Wtime();
  // End time

  cout << "The sum is " << sum << " and the time is " << finish - start << endl;
  return 0;
}
