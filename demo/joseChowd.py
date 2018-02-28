import json
import sys
sys.path.append("..");

import ast
import booleantools as bt
from booleantools import BooleanFunction
from itertools import permutations

permset = list(permutations([0,1,2,3,4]))

def getFuncList():
    with open("in/all", "r") as filein:
        out = [[bt.perms_orbit_polynomial(permset, BooleanFunction(ast.literal_eval(f), 5))] for f in filein]

    return out

def main():
    funcs = getFuncList()
    jose_chowder_set = bt.perms_orbit_polynomial(permset, BooleanFunction([[0],[1],[2],[3],[2,3]], 5))
    out_dict = {}

    for f_class in funcs:
        haus_dist = bt.hausdorff_distance_sets(jose_chowder_set, func)
        if out_dict[haus_dict] is not None:
            out_dict[haus_dict].append(str(f_class[0]))
        else:
            out_dict[haus_dict] = [f_class[0]]

    with open("out/jose-dist.json", "w") as fileout:
        fileout.write(json.dumps(out_dict))

main()
