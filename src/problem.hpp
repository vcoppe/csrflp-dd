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
        // problem constants
        double root_val;
        bitset<N> root_bits, terminal_bits;
        int n, *l, **c;

        // constraints
        int dep[N], prev[N], next[N];
        bitset<N> pred[N];

        // pre-computed lists for rough upper bound
        vector<tuple<int,int,int> > edges;
        vector<pair<int,int> > lengths;

        Problem(string filename);
        ~Problem();

        void add_constraints(string filename);

        // handle constraints
        bool feasible(shared_ptr<State> &parent, int var, int val);

        // dynamic programming model
        int successor(shared_ptr<State> &parent, int var, int val, shared_ptr<State> &child);

        // compute an upper bound based on edges and cuts given cut values and a set of vertices
        int combined_ub(int *state_cuts, bitset<N> &free, int n_free);
    private:
        bool read_file(string filename);
        double root_value();
};

#endif // PROBLEM_HPP
