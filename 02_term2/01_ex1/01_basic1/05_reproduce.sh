#!/bin/bash

mkdir -p out
gcc 01* -lm && ./a.out > out/01.txt
gcc 02* -lm && ./a.out > out/02.txt
gcc 03* -lm && ./a.out > out/03.txt
gcc 04* -lm && ./a.out > out/04.txt
