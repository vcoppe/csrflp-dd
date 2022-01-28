#include <algorithm>
#include <chrono>
#include <iostream>
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
int lowest_active_layer = 0, next_active_layer, best_lb = MIN_VALUE, best_ub = MAX_VALUE, best_sol[N], explored = 0, frontier_size = 0;
set<shared_ptr<FrontierNode>, FrontierNodeComparator> frontier[N+1];
set<shared_ptr<FrontierNode>, FrontierNodeComparator> ongoing;
unordered_map<bitset<N>,shared_ptr<FrontierNode> > mem[N+1];

// options
int max_width, max_time, n_threads, layer_cmp, step;
bool use_symmetry;

// synchronization
mutex global_mutex;

void log() {
    cout << "Explored " << explored
         << ", LB " << (problem->root_val - best_ub)
         << ", UB " << (problem->root_val - best_lb)
         << ", LAL " << lowest_active_layer
         << ", Frontier sz " << frontier_size << '\n';
}

bool must_stop() {
    if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start_time).count() > max_time) {
        cout << "Maximum time reached. Stopping thread.\n";
        return true;
    }
    return false;
}

int estimate(shared_ptr<Node> &node, int node_layer, shared_ptr<FrontierNode> &suffix_node) {
    shared_ptr<State> &state = node->state;
    int n_free = problem->n - node_layer;
    if (use_symmetry && n_free == lowest_active_layer) {
        unordered_map<bitset<N>,shared_ptr<FrontierNode> >::iterator mem_it = mem[lowest_active_layer].find(problem->root_bits & ~(state->free));
        if (mem_it != mem[lowest_active_layer].end()) {
            suffix_node = mem_it->second;
            return mem_it->second->lp;
        } else return MIN_VALUE;
    }
    return problem->combined_ub(state->cuts, state->free, n_free);
}

void clear(unordered_map<bitset<N>,shared_ptr<Node> > *graph) {
    for (int i=0; i<=problem->n; i++) graph[i].clear();
}

void develop(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root, int current_lb, int current_ub, shared_ptr<Node> &best_node, shared_ptr<FrontierNode> &best_suffix) {
    shared_ptr<Node> root_node = make_shared<Node>();
    root_node->state = root->state;
    root_node->lp = root->lp;
    root_node->ub = current_ub;

    graph[lowest_active_layer][root_node->state->free] = root_node;

    int node_idx;
    vector<shared_ptr<Node> > layer_nodes;

    for (int current_layer=lowest_active_layer; current_layer<problem->n; current_layer++) {
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
            if (current_layer == next_active_layer) parent->ub = MIN_VALUE; // reset to compute local bound
        }

        bool must_squash = current_layer+1 > next_active_layer && graph[current_layer+1].size() > max_width;

        if (must_squash) {
            if (layer_nodes.size() < graph[current_layer+1].size()) layer_nodes.resize(graph[current_layer+1].size());
            node_idx = 0;
        }

        for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[current_layer+1].begin(); it!=graph[current_layer+1].end(); ++it) {
            shared_ptr<Node> &node = it->second;
            shared_ptr<FrontierNode> suffix_node;
            int est = estimate(node, current_layer+1, suffix_node);
            node->ub = min(node->lp + est, current_ub);

            if (must_squash) layer_nodes[node_idx++] = node;
            if (suffix_node && node->ub > current_lb) {
                current_lb = node->ub;
                best_node = node;
                best_suffix = suffix_node;
            }

            if (current_layer == next_active_layer) {
                shared_ptr<Arc> arc = node->best_arc;
                while (arc != nullptr) { // quick local bound
                    arc->to->ub = max(arc->to->ub, arc->to->lp + arc->cost + est);
                    arc = arc->next_arc;
                }
            }
        }

        if (use_symmetry && problem->n-(current_layer+1) == lowest_active_layer) break;
        else if (must_squash) {
            int to_remove = graph[current_layer+1].size() - (lowest_active_layer ? max_width : 1000);
            if (to_remove > 0) {
                if (layer_cmp == 0) sort(layer_nodes.begin(), layer_nodes.begin()+node_idx, LayerNodeLPComparator());
                else sort(layer_nodes.begin(), layer_nodes.begin()+node_idx, LayerNodeUBComparator());
                for (node_idx=0; node_idx<to_remove; node_idx++) graph[current_layer+1].erase(layer_nodes[node_idx]->state->free);
            }
        }
    }

    if (!use_symmetry) {
        unordered_map<bitset<N>,shared_ptr<Node> >::iterator it = graph[problem->n].find(problem->terminal_bits);
        if (it != graph[problem->n].end()) best_node = it->second;
    }
}

void update_best(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root, shared_ptr<Node> &node) {
    if (node->lp > best_lb) {
        best_lb = node->lp;
        for (int layer=0; layer<lowest_active_layer; layer++) best_sol[layer] = root->sol[layer];
        bitset<N> cur;
        for (int layer=problem->n; layer>lowest_active_layer; layer--) {
            shared_ptr<Arc> &arc = graph[layer][cur]->best_arc;
            best_sol[layer-1] = arc->value;
            cur = arc->to->state->free;
        }
        log();
    }
}

void update_best_sym(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root, shared_ptr<Node> &prefix_node, shared_ptr<FrontierNode> &suffix_node) {
    if (prefix_node->lp + suffix_node->lp > best_lb) {
        best_lb = prefix_node->lp + suffix_node->lp;
        for (int layer=0; layer<lowest_active_layer; layer++) best_sol[layer] = root->sol[layer];
        bitset<N> cur = prefix_node->state->free;
        for (int layer=problem->n-lowest_active_layer; layer>lowest_active_layer; layer--) {
            shared_ptr<Arc> &arc = graph[layer][cur]->best_arc;
            best_sol[layer-1] = arc->value;
            cur = arc->to->state->free;
        }
        for (int layer=problem->n-lowest_active_layer; layer<problem->n; layer++) best_sol[layer] = suffix_node->sol[problem->n-layer-1];
        log();
    }
}

void enqueue_cutset(unordered_map<bitset<N>,shared_ptr<Node> > *graph, shared_ptr<FrontierNode> &root) {
    for (unordered_map<bitset<N>,shared_ptr<Node> >::iterator it=graph[next_active_layer].begin(); it!=graph[next_active_layer].end(); ++it) {
        shared_ptr<Node> &node = it->second;
        if (node->ub > best_lb) {
            // check if the state is already in frontier
            unordered_map<bitset<N>,shared_ptr<FrontierNode> >::iterator mem_it = mem[next_active_layer].find(node->state->free);
            if (mem_it != mem[next_active_layer].end()) { // state already there
                if (node->lp > mem_it->second->lp) { // if improves, update ub and reinsert in frontier
                    shared_ptr<FrontierNode> frontier_node = mem_it->second;
                    frontier[next_active_layer].erase(frontier[next_active_layer].find(frontier_node));
                    frontier_node->lp = node->lp;
                    frontier_node->ub = node->ub;

                    for (int layer=0; layer<lowest_active_layer; layer++) frontier_node->sol[layer] = root->sol[layer];
                    bitset<N> cur = node->state->free;
                    for (int layer=next_active_layer; layer>lowest_active_layer; layer--) {
                        shared_ptr<Arc> &arc = graph[layer][cur]->best_arc;
                        frontier_node->sol[layer-1] = arc->value;
                        cur = arc->to->state->free;
                    }

                    frontier[next_active_layer].insert(frontier_node);
                }
            } else {
                shared_ptr<FrontierNode> frontier_node = make_shared<FrontierNode>();
                frontier_node->state = node->state;
                frontier_node->lp = node->lp;
                frontier_node->ub = node->ub;
                frontier_node->layer = next_active_layer;

                for (int layer=0; layer<lowest_active_layer; layer++) frontier_node->sol[layer] = root->sol[layer];
                bitset<N> cur = node->state->free;
                for (int layer=next_active_layer; layer>lowest_active_layer; layer--) {
                    shared_ptr<Arc> &arc = graph[layer][cur]->best_arc;
                    frontier_node->sol[layer-1] = arc->value;
                    cur = arc->to->state->free;
                }

                frontier[next_active_layer].insert(frontier_node);
                mem[next_active_layer][node->state->free] = frontier_node;
                frontier_size++;
            }
        }
    }
}

void finalize_middle() {
    while (!frontier[(problem->n+1)/2].empty()) {
        set<shared_ptr<FrontierNode>, FrontierNodeComparator>::iterator it = --frontier[(problem->n+1)/2].end();
        shared_ptr<FrontierNode> prefix_node = *it;
        frontier[(problem->n+1)/2].erase(it);
        frontier_size--;

        unordered_map<bitset<N>,shared_ptr<FrontierNode> >::iterator mem_it = mem[problem->n/2].find(problem->root_bits & ~(prefix_node->state->free));
        if (mem_it != mem[problem->n/2].end()) { // find perfect match
            shared_ptr<FrontierNode> &suffix_node = mem_it->second;
            if (prefix_node->lp + suffix_node->lp > best_lb) {
                best_lb = prefix_node->lp + suffix_node->lp;
                for (int layer=0; layer<(problem->n+1)/2; layer++) best_sol[layer] = prefix_node->sol[layer];
                for (int layer=problem->n/2; layer<problem->n; layer++) best_sol[layer] = suffix_node->sol[problem->n-layer-1];
                log();
            }
        }
    }
}

void task() {
    unordered_map<bitset<N>,shared_ptr<Node> > graph[N+1];
    for (int i=0; i<=problem->n; i++) graph[i] = unordered_map<bitset<N>,shared_ptr<Node> >();

    while (true) {
        global_mutex.lock();

        if (must_stop() || lowest_active_layer == (problem->n+1)/2) { // time cutoff or end of algorithm
            global_mutex.unlock();
            break;
        } else if (frontier[lowest_active_layer].empty()) { // no more nodes in this layer
            if (ongoing.empty()) { // current layer fully developed
                if (next_active_layer == (problem->n+1)/2) { // layers ready to meet in the middle
                    finalize_middle();
                    mem[lowest_active_layer].clear();
                    mem[next_active_layer].clear();
                    lowest_active_layer = (problem->n+1)/2;
                } else {
                    mem[lowest_active_layer].clear();
                    lowest_active_layer = next_active_layer;
                    if (lowest_active_layer == problem->n/2 && problem->n%2 == 1) next_active_layer = (problem->n+1)/2;
                    else next_active_layer = min(lowest_active_layer+step, problem->n/2);
                }
                log();
            }
            global_mutex.unlock();
        } else {
            // get node to explore
            set<shared_ptr<FrontierNode>, FrontierNodeComparator>::iterator it = --frontier[lowest_active_layer].end();
            shared_ptr<FrontierNode> node = *it;
            frontier[lowest_active_layer].erase(it);
            frontier_size--;

            // update bounds
            int new_ub = node->ub;
            if (!ongoing.empty()) new_ub = max(new_ub, (*(--ongoing.end()))->ub);
            if (!frontier[next_active_layer].empty()) new_ub = max(new_ub, (*(--frontier[next_active_layer].end()))->ub);
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

            shared_ptr<Node> best_node;
            shared_ptr<FrontierNode> best_suffix;

            clear(graph);
            develop(graph, node, current_lb, current_ub, best_node, best_suffix);

            global_mutex.lock();

            if (!use_symmetry && best_node) update_best(graph, node, best_node);
            else if (best_node && best_suffix) update_best_sym(graph, node, best_node, best_suffix);
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

    frontier[0].insert(root);
    mem[0][root->state->free] = root;

    frontier_size++;

    next_active_layer = min(lowest_active_layer+step, problem->n/2);

    vector<thread> threads;
    for (int i=0; i<n_threads; i++) threads.emplace_back(task);
    for (int i=0; i<n_threads; i++) threads[i].join(); // wait until they finish
}

int main(int argc, char const *argv[]) {
    if (argc == 1 || (argc == 2 && strcmp(argv[1],"--help") == 0)) {
        cout << "Usage: ./PROGRAM filename [--width width]           // maximum width of the DDs\n"
             << "                          [--time maxtime]          // maximum time (seconds) for the algorithm\n"
             << "                          [--threads threads]       // number of threads used\n"
             << "                          [--symmetry-off]          // turn off symmetry-breaking\n"
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
    use_symmetry = 1;
    layer_cmp = 1;
    step = 1;

    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i],"--width") == 0) max_width = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--time") == 0) max_time = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--threads") == 0) n_threads = stoi(argv[i+1]);
        else if (strcmp(argv[i],"--symmetry-off") == 0) use_symmetry = 0;
        else if (strcmp(argv[i],"--minlp") == 0) layer_cmp = 0;
        else if (strcmp(argv[i],"--constraints") == 0) problem->add_constraints(argv[i+1]);
    }

    cout << "Starting search with width=" << max_width
         << ", time=" << max_time
         << ", n_threads=" << n_threads
         << ", exploit symmetry=" << (use_symmetry ? "on" : "off")
         << ", node comparator=" << (layer_cmp ? "ub" : "lp")
         << ", step=" << step << ".\n";

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
