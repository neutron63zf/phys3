#!/bin/bash

mkdir -p out


gcc 01* -lm -lblas -g -oa
./a

gcc 02* -lm -lblas -g -oa
./a