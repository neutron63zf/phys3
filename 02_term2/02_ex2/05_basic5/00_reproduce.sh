#!/bin/bash

gcc 01* -lm

# run $L $MS
function run() {
  L=$1
  MS=$2
  ignore=$3
  T=$4
  ./a.out $L $MS $ignore $T
}

L=2
MS=2
ignore=0
T=0.5

run $L $MS $ignore $T