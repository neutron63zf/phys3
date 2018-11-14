#!/bin/bash

## time系コマンドは標準出力に現れないので、時間だけはコピペが必要

mkdir -p out

# v -> height
v=10
# width -> width of center wall
width=0.1

# n -> division
n=100
./01* $n $v $width > out/01.txt

gcc 02* -lblas -lm

# runsim n offset mode
function runsim() {
  n=$1
  offset=$2
  mode=$3

  ./a.out $n $v $width $offset $mode
}

# runset n offset num1 num2
function runset() {
  n=$1
  offset=$2
  num1=$3
  num2=$4

  runsim $n $offset 0 > out/$num1.txt
  runsim $n $offset 1 > out/$num2.txt
}

# runtime n offset
function runtime() {
  n=$1
  offset=$2

  time runsim $n $offset 2
}

# offset -> potential offset

# mode=0 -> print eigenstate/value
# mode=1 -> print convergence
# mode=2 -> silent

n=25
offset=-1270
runtime $n $offset

n=50
offset=-5100
runtime $n $offset

n=100
offset=-21000
runtime $n $offset
runset $n $offset 02 03

n=150
offset=-45100
runtime $n $offset

n=200
offset=-81000
runtime $n $offset

n=250
offset=-125100
runtime $n $offset

n=300
offset=-181000
runtime $n $offset

n=350
offset=-245100
runtime $n $offset

n=400
offset=-321000
runtime $n $offset

n=500
offset=-510000
runtime $n $offset
