#!/bin/bash

# script to run the test cases and compare to reference cases

rm -rf *.mesh

for case in $(ls *.inp); do
  echo "====================================================="
  ../src/driver.exe $case > $case.out
done

for case in $(ls *.mesh); do
  echo "====================================================="
  echo "Diffing mesh for $case"
  diff $case.save $case
done


