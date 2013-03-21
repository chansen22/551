#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

int main() {
  int sum = 0;
  int tmp;
  clock_t t1, t2;
  vector<int> values;
  vector<int>::iterator it;
  while(cin >> tmp) {
    values.push_back(tmp);
  }


  // Start time
  t1 = clock();

  for(it = values.begin(); it != values.end(); ++it){
    sum = sum + *it;
  }

  // End time
  t2 = clock();
  float diff = (float)t2 - (float)t1;
  double seconds = diff / CLOCKS_PER_SEC;

  cout << "The sum is " << sum << " and the time is " << seconds << endl;
  return 0;
}
