CSCI 551 Numerical and Parallel Programming
-------------------------------------------

Compiling
---------

g++ -g -Wall -o serial_sum serial_sum.cpp
mpic++ -g -Wall -std=c++98 -o global_sum_send global_sum_send.cpp

Running
-------

./serial_sum < [file for input]
mpirun -np [number of cores] ./global_sum_send < [file for input]
