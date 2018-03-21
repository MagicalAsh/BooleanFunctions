from booleantools import BooleanFunction
import sys
import booleantools as bt
import itertools
import multiprocessing as mp
import json

# linear part + quadratic part
def analyze_polys(poly_gens, que, thread_no, n):
    for polys in poly_gens:
        for poly in polys:
            func = BooleanFunction(poly, n)
            if func.is_k_resilient(k=n-3) and not func.is_affine():
                que.put(poly)
    
    # Once all of the generators have executed
    que.put(thread_no)

def reduce_classes_dgen(class_list, new_f, n):
    """
    Produces a *MINIMALLY* reduced function list. This is by no means fully reduced.
    """
    for f in class_list:
        if new_f in f:
            return None

    class_list.append(get_class(new_f, n))
    return None

def reduce_classes(func_list):
    class_list = []
    for f in func_list:
        in_one = False
        for g in class_list:
            if f in g:
                in_one = True
        if in_one == False:
            class_list.append(f.get_orbit())
    return class_list

def get_class(f, n):
    perms = bt.Sym(n)
    return [apply_permutation(f, sigma) for sigma in perms]

def apply_permutation(poly, perm):
    def apply_perm_to_monomial(perm,monomial):
        out = []
        for var in monomial:
            out.append(perm[var])
        return out
       
    out = [apply_perm_to_monomial(perm, i) for i in poly]
    return out

def powerset(iterable):
    s = list(iterable)
    return list(itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)))

def func_generator(lin_part, non_lin_parts):
    lin = [[mon] for mon in lin_part]
    for nonlinear in non_lin_parts:
        yield lin + list(nonlinear)
            

def main():
    n = int(sys.argv[1])

    nonlin = powerset(itertools.combinations(list(range(n)), 2))
    nonlin.remove(())
    lin = powerset(list(range(n)))
    
    generator_list = [func_generator(linear, nonlin) for linear in lin]    
    size = len(generator_list)//16
    chunked = [generator_list[i:i+size] for i in range(0, len(generator_list), size)]         

    que = mp.Queue()
    threads = []
    for generators in chunked:
        thred = mp.Process(target=analyze_polys, args=(generators, que, len(threads), n))
        thred.start()
        threads.append(thred)
    
    deadCnt = 0
    f_out = []
    while deadCnt < len(chunked):
        f = que.get()
        if isinstance(f, list):
            reduce_classes_dgen(f_out, f, n)
        else: # it's an int
            deadCnt += 1
            threads[f].join()
    
    f_out = reduce_classes([BooleanFunction(f[0], n) for f in f_out])

    with open("out/classes_%dv.json" % n, "w") as outfile: 
        outfile.write(json.dumps({i: f_out[i][0].listform for i in range(len(f_out))}, indent=4))

if __name__ == "__main__":
    if(len(sys.argv) != 2):
        print("USAGE: python3 generate_classes.py <n>")
    else:
        main()

