# csrflp-dd
A Decision-Diagram-based approach to solve the Constrained Single Row Facility Layout Problem (cSRFLP).

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

The `data` folder contains benchmark instances for the Single-Row Facility Layout Problem (SRFLP), the Single Row Equidistant Facility Layout Problem (SREFLP) and the Minimum Linear Arrangement Problem (MinLA). The latter ones are special cases of the SRFLP. It also contains constraint files to transform any instance in an instance of the cSRFLP.

The file format for the SRFLP, SREFLP and MinLA instances is the following:
```
n
l_1  l_2  ... l_n
c_11 c_12 ... c_1n
c_21 c_22 ... c_2n
...  ...  ... ...
c_n1 c_n2 ... c_nn
```
Values can also be separated with commas.
The file format for the constraints is:

```
p o r
a_1 b_1 // position(a_1) = b_1
a_2 b_2
...
a_p b_p
c_1 d_1 // c_1 \in predecessors(d_1)
c_2 d_2
...
c_o d_o
e_1 f_1 // e_1 = previous(f_1)
e_2 f_2
...
e_r f_r
```
where `p`, `o` and `r` respectively denote the number of *positioning*, *ordering* and *relation* constraints.
The following lines contain unique constraints on departments and positions, both 0-indexed.
