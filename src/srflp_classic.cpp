#include <algorithm>
#include <chrono>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include <thread>
#include <unordered_map>

#include "structs.hpp"
#include "problem.hpp"

unique_ptr<Problem> problem;

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
         << ", LB " << (problem->root_val - best_ub)
         << ", UB " << (problem->root_val - best_lb)
         << ", Frontier sz " << frontier.size() << '\n';
}

bool must_stop() {
    if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start_time).count() > max_time) {
        cout << "Maximum time reached. Stopping thread.\n";
        return true;
    }
    return false;
}

int estimate(shared_ptr<Node> &node, int node_layer) {
    shared_ptr<State> &state = node->state;
    int n_free = problem->n - node_layer;
    return problem->combined_ub(state->cuts, state->free, n_free);
}

void clear(unordered_map<bitset<N>,shared_ptr<Node> > *graph) {
    for (int i=0; i<=problem->n; i++) graph[i].clear();
}

void develop(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root, int current_lb, int current_ub) {
    shared_ptr<Node> root_node = make_shared<Node>();
    root_node->state = root->state;
    root_node->lp = root->lp;
    root_node->ub = current_ub;

    int first_layer = root->layer, node_idx;
    vector<shared_ptr<Node> > layer_nodes;

    graph[first_layer][root_node->state->free] = root_node;

    for (int current_layer=first_layer; current_layer<problem->n; current_layer++) {
        for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[current_layer].begin(); it!=graph[current_layer].end(); ++it) {
            shared_ptr<Node> &parent = it->second;
            if (parent->ub > current_lb) {
                for (int val=0; val<problem->n; val++) if (parent->state->free[val] && problem->feasible(parent->state, current_layer, val)) {
                    shared_ptr<State> child = nullptr;
                    int cost = problem->successor(parent->state, current_layer, val, child);

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
            int to_remove = graph[current_layer+1].size() - max_width;
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
    unordered_map<bitset<N>,shared_ptr<Node> >::iterator it = graph[problem->n].find(problem->terminal_bits);
    if (it != graph[problem->n].end() && it->second->lp > best_lb) {
        best_lb = it->second->lp;
        for (int layer=0; layer<first_layer; layer++) best_sol[layer] = root->sol[layer];
        bitset<N> cur;
        for (int layer=problem->n; layer>first_layer; layer--) {
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
    for (int i=0; i<=problem->n; i++) graph[i] = unordered_map<bitset<N>,shared_ptr<Node> >();

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
    root->state->free = problem->root_bits;

    frontier.insert(root);
    mem[root->state->free] = root;

    vector<thread> threads;
    for (int i=0; i<n_threads; i++) threads.emplace_back(task);
    for (int i=0; i<n_threads; i++) threads[i].join(); // wait until they finish
}

int main(int argc, char const *argv[]) {
    if (argc == 1 || (argc == 2 && strcmp(argv[1],"--help") == 0)) {
        cout << "Usage: ./PROGRAM filename [--width width]           // maximum width of the DDs\n"
             << "                          [--time maxtime]          // maximum time (seconds) for the algorithm\n"
             << "                          [--threads threads]       // number of threads used\n"
             << "                          [--minlp]                 // use MinLP for comparison before restriction\n"
             << "                          [--constraints filename]  // add constraints to the problem\n";
        return 0;
    }

    problem = make_unique<Problem>(argv[1]);

    if (problem->n == -1) return 0; // error

    cout.precision(10);

    // default args
    max_width = problem->n;
    max_time = MAX_VALUE;
    n_threads = thread::hardware_concurrency();
    layer_cmp = 1;

    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i],"--width") == 0) max_width = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--time") == 0) max_time = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--threads") == 0) n_threads = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--minlp") == 0) layer_cmp = 0;
        else if (strcmp(argv[i],"--constraints") == 0) problem->add_constraints(argv[i+1]);
    }

    cout << "Starting search with width=" << max_width
         << ", time=" << max_time
         << ", n_threads=" << n_threads
         << ", node comparator=" << (layer_cmp ? "ub" : "lp") << ".\n";

    start_time = chrono::high_resolution_clock::now();
    solve();
    end_time = chrono::high_resolution_clock::now();

    cout << "Best value found: " << (problem->root_val - best_lb) << '\n';
    cout << "Layout          : ";
    for (int i=0; i<problem->n; i++) cout << best_sol[i] << " ";
    cout << "\nTime elapsed    : " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << "ms\n";
    cout << "Nodes explored  : " << explored << '\n';

    return 0;
}
