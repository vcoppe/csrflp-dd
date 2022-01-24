#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <tuple>

#include "structs.hpp"

class Problem {
    public:
        double root_val;
        bitset<N> root_bits, terminal_bits;
        int n, *l, **c;
        vector<tuple<int,int,int> > edges;
        vector<pair<int,int> > lengths;

        Problem(string filename);
        ~Problem();

        // dynamic programming model
        int successor(shared_ptr<State> &parent, int i, shared_ptr<State> &child);

        // compute an upper bound based on edges and cuts given cut values and a set of vertices
        int combined_ub(int *state_cuts, bitset<N> &free, int n_free);
    private:
        bool read_file(string filename);
        double root_value();
};

#endif // PROBLEM_HPP
