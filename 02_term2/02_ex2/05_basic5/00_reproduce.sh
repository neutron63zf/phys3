#!/bin/bash

mkdir -p out
gcc 01* -lm -g

# run $L $MS $ignore $T
function run() {
  L=$1
  MS=$2
  ignore=$3
  T=$4
  ./a.out $L $MS $ignore $T 0
}

# trans $L $MS $ignore
function trans() {
  L=$1
  MS=$2
  ignore=$3
  ./a.out $L $MS $ignore -1 1
}

L=8
ignore=10000
MS=$(($ignore * 10))
T=2

ignore=0
#run $L $MS $ignore $T > ./out/01.txt

ignore=10000
#run $L $MS $ignore $T > ./out/02.txt

L=12
#run $L $MS $ignore $T > ./out/03.txt

L=16
#run $L $MS $ignore $T > ./out/04.txt


T=3
L=8
#run $L $MS $ignore $T > ./out/05.txt

L=12
#run $L $MS $ignore $T > ./out/06.txt

L=16
#run $L $MS $ignore $T > ./out/07.txt

ignore=2500
MS=$(($ignore * 10))

L=8
trans $L $MS $ignore > ./out/08.txt &

L=12
trans $L $MS $ignore > ./out/09.txt &

L=16
trans $L $MS $ignore > ./out/10.txt &