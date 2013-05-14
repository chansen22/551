#!/bin/bash

for cores in 1 2 5 10 20 30
do
  for run in 1 2 3 4 5
  do 
    out=`./matrix_sum 8000 $cores`
    echo "Run $run with $cores cores the time was:"
    echo $out
  done
done
