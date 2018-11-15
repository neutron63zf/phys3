#!/bin/bash

mkdir -p out
gcc 01* -lm

# run $L $MS
function run() {
  L=$1
  MS=$2
  ignore=$3
  T=$4
  ./a.out $L $MS $ignore $T
}

L=8
ignore=1000
MS=$(($ignore * 10))
T=2

ignore=0
run $L $MS $ignore $T > ./out/01.txt

ignore=1000
run $L $MS $ignore $T > ./out/02.txt

L=12
run $L $MS $ignore $T > ./out/03.txt

L=16
run $L $MS $ignore $T > ./out/04.txt