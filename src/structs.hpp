#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include <bitset>

using namespace std;

// constants
#define MIN_VALUE -2147483648
#define MAX_VALUE 2147483647

// problem data
#define N 100 // n <= N

// mdd data structures
typedef struct {
    bitset<N> free;
    int cuts[N];
} State;

struct Node;

typedef struct Arc {
    shared_ptr<Node> to;
    shared_ptr<Arc> next_arc;
    int cost;
    int value;
} Arc;

typedef struct Node {
    shared_ptr<State> state;
    int lp;
    int ub;
    shared_ptr<Arc> best_arc;
} Node;

typedef struct {
    shared_ptr<State> state;
    int lp;
    int ub;
    int layer;
    int sol[N];
} FrontierNode;

// comparators
struct FrontierNodeComparator {
    inline bool operator() (const shared_ptr<FrontierNode> &lhs, const shared_ptr<FrontierNode> &rhs) const {
        if (lhs->ub == rhs->ub) {
            if (lhs->lp == rhs->lp) {
                for (int i=0; i<N; i++) {
                    if (lhs->state->free[i] ^ rhs->state->free[i]) return lhs->state->free[i];
                }
                return false;
            }
            return lhs->lp < rhs->lp;
        }
        return lhs->ub < rhs->ub;
    }
};

struct LayerNodeLPComparator {
    inline bool operator() (const shared_ptr<Node> &lhs, const shared_ptr<Node> &rhs) {
        if (lhs->lp == rhs->lp) return lhs->ub < rhs->ub;
        return lhs->lp < rhs->lp;
    }
};

struct LayerNodeUBComparator {
    inline bool operator() (const shared_ptr<Node> &lhs, const shared_ptr<Node> &rhs) {
        if (lhs->ub == rhs->ub) return lhs->lp < rhs->lp;
        return lhs->ub < rhs->ub;
    }
};

#endif // STRUCTS_HPP
