# srflp-dd
A Decision-Diagram-based approach to solve the Single Row Facility Layout Problem (SRFLP).

* `srflp_classic.cpp` performs a best-first branch-and-bound.
* `srflp_boosted.cpp` performs a breadth-first branch-and-bound and uses a symmetry-breaking technique (turned off by default because not valid with constraints).

## Build

If you have `g++` installed, just run `./compile.sh src/srflp_classic.cpp` and `./compile.sh src/srflp_boosted.cpp`.

## Run

`srflp_classic.cpp`
```
Usage: ./srflp_classic filename [--width width]          // maximum width of the DDs
                                [--time maxtime]         // maximum time (seconds) for the algorithm
                                [--threads threads]      // number of threads used
                                [--minlp]                // use MinLP for comparison before restriction
                                [--constraints filename] // add constraints to the problem
```

`srflp_boosted.cpp`
```
Usage: ./srflp_boosted filename [--width width]          // maximum width of the DDs
                                [--time maxtime]         // maximum time (seconds) for the algorithm
                                [--threads threads]      // number of threads used
                                [--symmetry-on]          // turn on symmetry-breaking
                                [--minlp]                // use MinLP for comparison before restriction
                                [--constraints filename] // add constraints to the problem
```

## Data

The `data` folder contains benchmark instances for the SRFLP, the Single Row Equidistant Facility Layout Problem (SREFLP) and the Minimum Linear Arrangement Problem (MinLA). The latter ones are special cases of the SRFLP.
