#!/bin/bash

g++ $1 -o $(basename $1 .cpp) -O3 -std=c++11 -pthread
