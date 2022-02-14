#!/usr/bin/python

import sys, getopt, random, time
from ortools.sat.python import cp_model
from pathlib import Path

random.seed(42)

n = 0
p = 0
o = 0
r = 0
k = 1

pc = None
oc = None
rc = None

last_model = None

solver = cp_model.CpSolver()

def init():
    global last_model, pc, oc, rc

    last_model = cp_model.CpModel()
    x = [last_model.NewIntVar(0, n-1, 'x%i' % i) for i in range(n)]
    last_model.AddAllDifferent(x)

    pc = {}
    oc = set()
    rc = {}

def model_with_additional_positional_constraint(dep, pos):
    model = cp_model.CpModel()
    model.CopyFrom(last_model)
    model.Add(model.GetIntVarFromProtoIndex(dep) == pos)
    return model

def model_with_additional_ordering_constraint(dep1, dep2):
    model = cp_model.CpModel()
    model.CopyFrom(last_model)
    model.Add(model.GetIntVarFromProtoIndex(dep1) < model.GetIntVarFromProtoIndex(dep2))
    return model

def model_with_additional_relational_constraint(dep1, dep2):
    model = cp_model.CpModel()
    model.CopyFrom(last_model)
    model.Add(model.GetIntVarFromProtoIndex(dep1)+1 == model.GetIntVarFromProtoIndex(dep2))
    return model

def feasible(model, print_sol=False):
    global last_model

    status = solver.Solve(model)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        last_model = model
        if print_sol:
            for i in range(n):
                print(solver.Value(model.GetIntVarFromProtoIndex(i)))
        return True
    else:
        return False

def generate():
    while len(pc) < p:
        dep = random.randint(0, n-1)
        pos = random.randint(0, n-1)

        if dep not in pc.keys() and pos not in pc.values():
            if feasible(model_with_additional_positional_constraint(dep, pos)):
                pc[dep] = pos

    while len(oc) < o:
        dep1 = random.randint(0, n-1)
        dep2 = random.randint(0, n-1)

        if dep1 != dep2 and (dep1,dep2) not in oc:
            if feasible(model_with_additional_ordering_constraint(dep1,dep2)):
                oc.add((dep1,dep2))

    while len(rc) < r:
        dep1 = random.randint(0, n-1)
        dep2 = random.randint(0, n-1)

        if dep1 not in rc.keys() and dep2 not in rc.values():
            if feasible(model_with_additional_relational_constraint(dep1,dep2)):
                rc[dep1] = dep2

    #feasible(last_model,True)

def write(i):
    Path("./constraints/%d" % n).mkdir(parents=True, exist_ok=True)
    f = open("./constraints/%d/srflp_%d_%d_%d_%d" % (n, p, o, r, i), "w")
    f.write("%d %d %d\n" % (p, o, r))
    for dep, pos in pc.items():
        f.write("%d %d\n" % (dep, pos))
    for (dep1,dep2) in oc:
        f.write("%d %d\n" % (dep1, dep2))
    for dep1, dep2 in rc.items():
        f.write("%d %d\n" % (dep1, dep2))
    f.close()

def main(argv):
    global n, p, o, r, k

    input_format = """usage:
    srflp_constraints.py
    -n <number of departments>
    -p <number of positional constraints>
    -o <number of ordering constraints>
    -r <number of relational constraints>
    -k <number of instances to generate>"""

    try:
        opts, args = getopt.getopt(argv,"hn:p:o:r:k:")
        for opt, arg in opts:
            match opt:
                case '-h':
                    print(input_format)
                    sys.exit()
                case '-n': n = int(arg)
                case '-p': p = int(arg)
                case '-o': o = int(arg)
                case '-r': r = int(arg)
                case '-k': k = int(arg)
    except getopt.GetoptError:
        print(input_format)
        sys.exit(2)

    if n < 1:
        print('error: n must be strictly positive')
        sys.exit(2)
    if p > n:
        print('error: p must be smaller or equal to n')
        sys.exit(2)
    if o > (n * (n-1)) / 2:
        print('error: o must be smaller or equal to the total number of possible pairs')
        sys.exit(2)
    if r >= n:
        print('error: r must be smaller than n')
        sys.exit(2)
    if k < 1:
        print('error: k must be striclty positive')
        sys.exit(2)

    for i in range(k):
        init()
        generate()
        write(i)

if __name__ == "__main__":
    main(sys.argv[1:])
