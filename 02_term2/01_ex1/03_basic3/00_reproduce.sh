#!/bin/bash

mkdir -p out

# division
n=100
# height
v=10
width=0.1
offset=-21000
# mode=1 -> print convergence
# mode=0 -> print eigenstate/value
mode=0

#./01* $n $v $width #> out/01.txt
gcc 02* -lblas -lm && time ./a.out $n $v $width $offset $mode #> out/02.txt
