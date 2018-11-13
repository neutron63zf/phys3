#!/bin/bash

mkdir -p out

# division
n=200
# height
v=0
width=0.1
offset=-160000

./01* $n $v $width #> out/01.txt
#gcc 02* -lblas -lm && time ./a.out $n $v $width $offset #> out/02.txt
