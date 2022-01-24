#include "problem.hpp"

Problem::Problem(string filename) {
    if (!read_file(filename)) {
        n = -1;
        return;
    }

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

Problem::~Problem() {
    if (n == -1) return;

    delete [] l;
    for (int i=0; i<n; i++) delete [] c[i];
    delete [] c;
}

bool Problem::read_file(string filename) {
    string line, number;
    smatch match;
    ifstream file(filename);
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
                        return false;
                    }
                    l = new int[n];
                    c = new int*[n];
                    for (int i=0; i<n; i++) c[i] = new int[n];
                } else if (cnt > 0 && cnt <= n) {
                    l[cnt-1] = stoi(number);
                    if (regex_search(filename, clearance)) l[cnt-1] += 10; // instances with clearance requirements
                } else if (cnt > n && cnt <= n + n*n) {
                    c[(cnt-n-1)/n][(cnt-n-1)%n] = stoi(number);
                }
                line = match.suffix().str();
                cnt++;
            }
        }
        file.close();
        return true;
    } else {
        cout << "Error opening file.\n";
        return false;
    }
}

double Problem::root_value() {
    double offset = 0;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i != j) {
        offset += 0.5 * c[i][j] * l[i];
    }
    return offset;
}

int Problem::successor(shared_ptr<State> &parent, int i, shared_ptr<State> &child) {
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

int Problem::combined_ub(int *state_cuts, bitset<N> &free, int n_free) {
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
