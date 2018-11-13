#!/bin/bash

mkdir -p out
gcc 01* -lm && ./a.out > out/01.txt
