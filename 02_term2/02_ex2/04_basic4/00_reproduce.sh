#!/bin/bash

mkdir -p out

gcc 01* -lm

p_c=0.1
lambda=5
M=10000
D_start=1
D_end=11
D_step_num=100

# run $p_c $lambda $M $D_start $D_end $D_step_num
function run() {
  p_c=$1
  lambda=$2
  M=$3
  D_start=$4
  D_end=$5
  D_step_num=$6

  ./a.out $p_c $lambda $M $D_start $D_end $D_step_num
}

run $p_c $lambda $M $D_start $D_end $D_step_num > ./out/01.txt