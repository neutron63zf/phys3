#!/bin/bash

mkdir -p out
gcc 02* -lm -g

# run $L $T
function run() {
  L=$1
  T=$2
  ./a.out $L $T
}

run 4 2.0000

echo ""
run 5 2.0000

# todo