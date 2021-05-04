#include <algorithm>
#include <bitset>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <regex>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>

using namespace std;

// constants
#define MIN_VALUE -2147483648
#define MAX_VALUE 2147483647

// problem data
#define N 100 // n <= N
int n, l[N], c[N][N];
vector<tuple<int,int,int> > edges;
vector<pair<int,int> > lengths;
bitset<N> root_bits, terminal_bits;
double root_val;

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
                for (int i=0; i<n; i++) {
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

// search variables
chrono::time_point<chrono::high_resolution_clock> start_time, end_time;
int best_lb = MIN_VALUE, best_ub = MAX_VALUE, best_sol[N], explored = 0;
set<shared_ptr<FrontierNode>, FrontierNodeComparator> frontier, ongoing;
unordered_map<bitset<N>,shared_ptr<FrontierNode> > mem;

// options
int max_width, max_time, layer_cmp, n_threads;

// synchronization
mutex global_mutex;

void log() {
    cout << "Explored " << explored
         << ", LB " << (root_val - best_ub)
         << ", UB " << (root_val - best_lb)
         << ", Frontier sz " << frontier.size() << '\n';
}

bool must_stop() {
    if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start_time).count() > max_time) {
        cout << "Maximum time reached. Stopping thread.\n";
        return true;
    }
    return false;
}

double root_value() {
    double offset = 0;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i != j) {
        offset += 0.5 * c[i][j] * l[i];
    }
    return offset;
}

// compute an upper bound based on edges and cuts given cut values and a set of vertices
int combined_ub(int *state_cuts, bitset<N> &free, int n_free) {
    vector<pair<double, int> > order;

    int edge_lb = 0, cut_lb = 0, cumul_l = 0, edge_idx = 0, length_idx = 0;
    for (int i=0; i<n_free; i++) {
        for (int k=0; k<n_free-i-1; k++) {
            while (!free[get<1>(edges[edge_idx])] || !free[get<2>(edges[edge_idx])]) edge_idx++;
            edge_lb += cumul_l * get<0>(edges[edge_idx]);
            edge_idx++;
        }

        while (!free[lengths[length_idx].second]) length_idx++;
        cumul_l += lengths[length_idx].first;

        order.emplace_back(- ((double) state_cuts[lengths[length_idx].second]) / lengths[length_idx].first, lengths[length_idx].second);

        length_idx++;
    }

    sort(order.begin(), order.end());

    cumul_l = 0;
    for (int i=0; i<order.size(); i++) {
        int j = order[i].second;
        cut_lb += cumul_l * state_cuts[j];
        cumul_l += l[j];
    }

    return - edge_lb - cut_lb;
}

int estimate(shared_ptr<Node> &node, int node_layer) {
    shared_ptr<State> &state = node->state;
    int n_free = n - node_layer;
    return combined_ub(state->cuts, state->free, n_free);
}

// dynamic programming model
int successor(shared_ptr<State> &parent, int i, shared_ptr<State> &child) {
    child = make_shared<State>();
    child->free |= parent->free;
    child->free[i] = false;
    child->cuts[i] = 0;

    int cut = 0;
    for (int j=0; j<n; j++) if (child->free[j]) {
        cut += parent->cuts[j];
        child->cuts[j] = parent->cuts[j] + c[i][j];
    }

    return - l[i] * cut;
}

void init() {
    for (int i=0; i<n; i++) root_bits[i] = true;

    // precompute lists for estimate
    for (int i=0; i<n; i++) {
        lengths.push_back(make_pair(l[i], i));
        for (int j=i+1; j<n; j++) {
            if (c[i][j] != c[j][i]) c[i][j] = c[j][i] = c[i][j] + c[j][i];
            edges.push_back(make_tuple(c[i][j], i, j));
        }
    }

    root_val = root_value();

    sort(lengths.begin(), lengths.end());
    sort(edges.begin(), edges.end(), greater<tuple<int,int,int> >());
}

void clear(unordered_map<bitset<N>,shared_ptr<Node> > *graph) {
    for (int i=0; i<=n; i++) graph[i].clear();
}

void develop(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root, int current_lb, int current_ub) {
    shared_ptr<Node> root_node = make_shared<Node>();
    root_node->state = root->state;
    root_node->lp = root->lp;
    root_node->ub = current_ub;

    int first_layer = root->layer, node_idx;
    vector<shared_ptr<Node> > layer_nodes;

    graph[first_layer][root_node->state->free] = root_node;

    for (int current_layer=first_layer; current_layer<n; current_layer++) {
        for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[current_layer].begin(); it!=graph[current_layer].end(); ++it) {
            shared_ptr<Node> &parent = it->second;
            if (parent->ub > current_lb) {
                for (int val=0; val<n; val++) if (parent->state->free[val]) {
                    shared_ptr<State> child = nullptr;
                    int cost = successor(parent->state, val, child);

                    shared_ptr<Arc> arc = make_shared<Arc>();
                    arc->to = parent;
                    arc->cost = cost;
                    arc->value = val;

                    unordered_map<bitset<N>,shared_ptr<Node> >::iterator node_it=graph[current_layer+1].find(child->free);
                    if (node_it != graph[current_layer+1].end()) {
                        shared_ptr<Node> &node = node_it->second;
                        if (parent->lp + cost > node->lp) {
                            node->lp = parent->lp + cost;
                            arc->next_arc = node->best_arc;
                            node->best_arc = arc;
                        } else {
                            arc->next_arc = node->best_arc->next_arc;
                            node->best_arc->next_arc = arc;
                        }
                    } else {
                        shared_ptr<Node> node = make_shared<Node>();
                        node->state = child;
                        node->best_arc = arc;
                        node->lp = parent->lp + cost;
                        node->ub = current_ub;

                        graph[current_layer+1][node->state->free] = node;
                    }
                }
            }
            if (current_layer == first_layer+1) parent->ub = MIN_VALUE; // reset to compute local bound
        }

        bool must_squash = current_layer > first_layer && graph[current_layer+1].size() > max_width;

        if (must_squash) {
            if (layer_nodes.size() < graph[current_layer+1].size()) layer_nodes.resize(graph[current_layer+1].size());
            node_idx = 0;
        }

        for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[current_layer+1].begin(); it!=graph[current_layer+1].end(); ++it) {
            shared_ptr<Node> &node = it->second;
            shared_ptr<FrontierNode> suffix_node;
            int est = estimate(node, current_layer+1);
            node->ub = min(node->lp + est, current_ub);

            if (must_squash) layer_nodes[node_idx++] = node;

            if (current_layer == first_layer+1) {
                shared_ptr<Arc> arc = node->best_arc;
                while (arc != nullptr) { // quick local bound
                    arc->to->ub = max(arc->to->ub, arc->to->lp + arc->cost + est);
                    arc = arc->next_arc;
                }
            }
        }

        if (must_squash) {
            int to_remove = graph[current_layer+1].size() - (root->layer ? max_width : 1000);
            if (to_remove > 0) {
                if (layer_cmp == 0) sort(layer_nodes.begin(), layer_nodes.begin()+node_idx, LayerNodeLPComparator());
                else sort(layer_nodes.begin(), layer_nodes.begin()+node_idx, LayerNodeUBComparator());
                for (node_idx=0; node_idx<to_remove; node_idx++) graph[current_layer+1].erase(layer_nodes[node_idx]->state->free);
            }
        }
    }
}

void update_best(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root) {
    int first_layer = root->layer;
    unordered_map<bitset<N>,shared_ptr<Node> >::iterator it = graph[n].find(terminal_bits);
    if (it != graph[n].end() && it->second->lp > best_lb) {
        best_lb = it->second->lp;
        for (int layer=0; layer<first_layer; layer++) best_sol[layer] = root->sol[layer];
        bitset<N> cur;
        for (int layer=n; layer>first_layer; layer--) {
            shared_ptr<Arc> &arc = graph[layer][cur]->best_arc;
            best_sol[layer-1] = arc->value;
            cur = arc->to->state->free;
        }
        log();
    }
}

void enqueue_cutset(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root) {
    int first_layer = root->layer;
    for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[first_layer+1].begin(); it!=graph[first_layer+1].end(); ++it) {
        shared_ptr<Node> &node = it->second;
        if (node->ub > best_lb) {
            // check if the state is already in frontier
            unordered_map<bitset<N>,shared_ptr<FrontierNode> >::iterator mem_it = mem.find(node->state->free);
            if (mem_it != mem.end()) { // state already there
                if (node->lp > mem_it->second->lp) { // if improves, update ub and reinsert in frontier
                    shared_ptr<FrontierNode> frontier_node = mem_it->second;
                    frontier.erase(frontier.find(frontier_node));
                    frontier_node->lp = node->lp;
                    frontier_node->ub = node->ub;

                    for (int layer=0; layer<first_layer; layer++) frontier_node->sol[layer] = root->sol[layer];
                    frontier_node->sol[first_layer] = graph[first_layer+1][node->state->free]->best_arc->value;

                    frontier.insert(frontier_node);
                }
            } else {
                shared_ptr<FrontierNode> frontier_node = make_shared<FrontierNode>();
                frontier_node->state = node->state;
                frontier_node->lp = node->lp;
                frontier_node->ub = node->ub;
                frontier_node->layer = first_layer+1;

                for (int layer=0; layer<first_layer; layer++) frontier_node->sol[layer] = root->sol[layer];
                frontier_node->sol[first_layer] = graph[first_layer+1][node->state->free]->best_arc->value;

                frontier.insert(frontier_node);
                mem[node->state->free] = frontier_node;
            }
        }
    }
}

void task() {
    unordered_map<bitset<N>,shared_ptr<Node> > graph[N+1];
    for (int i=0; i<=n; i++) graph[i] = unordered_map<bitset<N>,shared_ptr<Node> >();

    while (true) {
        global_mutex.lock();

        if (must_stop()) { // time cutoff
            global_mutex.unlock();
            break;
        } else if (frontier.empty()) {
            if (ongoing.empty()) { // end of algorithm
                global_mutex.unlock();
                break;
            } else { // waiting for other threads to finish or enqueue other nodes
                global_mutex.unlock();
            }
        } else {
            // get node to explore
            set<shared_ptr<FrontierNode>, FrontierNodeComparator>::iterator it = --frontier.end();
            shared_ptr<FrontierNode> node = *it;
            frontier.erase(it);
            mem.erase(node->state->free);

            if (ongoing.find(node) != ongoing.end()) { // node already being developed
                global_mutex.unlock();
                continue;
            }

            // update bounds
            int new_ub = node->ub;
            if (!ongoing.empty()) new_ub = max(new_ub, (*(--ongoing.end()))->ub);
            best_ub = min(best_ub, new_ub);

            if (best_lb >= best_ub) { // optimal value found
                global_mutex.unlock();
                break;
            }

            if (node->ub < best_lb) { // can be pruned
                global_mutex.unlock();
                continue;
            }

            // get current bounds
            int current_lb = best_lb, current_ub = best_ub;

            explored++; // increment nodes explored
            if (explored % 10000 == 0) log();

            ongoing.insert(node);

            global_mutex.unlock();

            clear(graph);
            develop(graph, node, current_lb, current_ub);

            global_mutex.lock();

            update_best(graph, node);
            enqueue_cutset(graph, node);

            ongoing.erase(ongoing.find(node));

            global_mutex.unlock();
        }
    }
}

void solve() {
    shared_ptr<FrontierNode> root = make_shared<FrontierNode>();
    root->state = make_shared<State>();
    root->ub = MAX_VALUE;
    root->state->free = root_bits;

    frontier.insert(root);
    mem[root->state->free] = root;

    vector<thread> threads;
    for (int i=0; i<n_threads; i++) threads.emplace_back(task);
    for (int i=0; i<n_threads; i++) threads[i].join(); // wait until they finish
}

int main(int argc, char const *argv[]) {
    if (argc == 1 || (argc == 2 && strcmp(argv[1],"--help") == 0)) {
        cout << "Usage: ./PROGRAM filename [--width width]     // maximum width of the DDs\n"
             << "                          [--time maxtime]    // maximum time (seconds) for the algorithm\n"
             << "                          [--threads threads] // number of threads used\n"
             << "                          [--minlp]           // use MinLP for comparison before restriction\n";
        return 0;
    }

    string line, number;
    smatch match;
    ifstream file(argv[1]);
    if (file.is_open()) {
        regex digit("(\\d+)"), non_digit("(\\D+)"), clearance("Cl");
        int cnt = 0;
        while (getline(file,line)) {
            line = regex_replace(line, non_digit, " ");
            while (regex_search(line, match, digit)) {
                number = match.str(0);
                if (cnt == 0) {
                    n = stoi(number);
                    if(n > N) {
                        cout << "Increase value of N in implementation.\n";
                        return 0;
                    }
                } else if (cnt > 0 && cnt <= n) {
                    l[cnt-1] = stoi(number);
                    if (regex_search(argv[1], clearance)) l[cnt-1] += 10; // instances with clearance requirements
                } else if (cnt > n && cnt <= n + n*n) {
                    c[(cnt-n-1)/n][(cnt-n-1)%n] = stoi(number);
                }
                line = match.suffix().str();
                cnt++;
            }
        }
        file.close();
    } else {
        cout << "Error opening file.\n";
        return 0;
    }

    cout.precision(10);

    // default args
    max_width = n;
    max_time = MAX_VALUE;
    n_threads = thread::hardware_concurrency();
    layer_cmp = 1;

    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i],"--width") == 0) max_width = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--time") == 0) max_time = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--threads") == 0) n_threads = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--minlp") == 0) layer_cmp = 0;
    }

    cout << "Starting search with width=" << max_width
         << ", time=" << max_time
         << ", n_threads=" << n_threads
         << ", node comparator=" << (layer_cmp ? "ub" : "lp") << ".\n";

    init();

    start_time = chrono::high_resolution_clock::now();
    solve();
    end_time = chrono::high_resolution_clock::now();

    cout << "Best value found: " << (root_val - best_lb) << '\n';
    cout << "Layout          : ";
    for (int i=0; i<n; i++) cout << best_sol[i] << " ";
    cout << "\nTime elapsed    : " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << "ms\n";
    cout << "Nodes explored  : " << explored << '\n';

    return 0;
}
