#!/bin/bash

g++ src/problem.cpp $1 -o $(basename $1 .cpp) -O3 -std=c++14 -pthread
