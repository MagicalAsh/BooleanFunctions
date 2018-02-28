from booleantools import BooleanFunction
import booleantools as bt
import itertools
import queue
import threading
import json

# linear part + quadratic part
def analyze_polys(polys, que):
    for poly in polys:
        if poly.is_k_resilient(k=2) and not poly.is_affine():
            que.put(poly)
    que.put(None)

def reduce_classes(func_list):
    class_list = []
    for f in func_list:
        in_one = False
        for g in class_list:
            if f in g:
                in_one = True
        if in_one == False:
            class_list.append(f.get_class())
    return class_list

def powerset(iterable):
    s = list(iterable)
    return list(itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)))

def main():
    nonlin = powerset(itertools.combinations([0,1,2,3,4], 2))
    nonlin.remove(())
    lin = powerset([0,1,2,3,4])
    f_list = []
    for nonlinear in nonlin:
        for linear in lin:
            func = []
            for i in linear:
                func.append([i])
            for i in nonlinear:
                func.append(list(i))
            f_list.append(BooleanFunction(func, 5))
    
    print(len(f_list))
    chunked = [f_list[i:i+500] for i in range(0, len(f_list), 500)]         
    que = queue.Queue()
    for lst in chunked:
        thred = threading.Thread(target=analyze_polys, args=(lst, que))
        thred.start()
    
    noneCnt = 0
    f_out = []
    while noneCnt < len(chunked):
        f = que.get()
        if f is not None:
            f_out.append(f)
        else:
            noneCnt += 1
    
    print("reducing classes: there are %d functions to reduce" % (len(f_out)))
    classes = reduce_classes(f_out)
    for f_class in classes:
        print(f_class[0].__repr__())
    
    with open("out/classes.json", "w") as outfile:
        outfile.write(json.dumps({i: classes[i][0].listform for i in range(len(classes))}, indent=4))

if __name__ == "__main__":
    main()

