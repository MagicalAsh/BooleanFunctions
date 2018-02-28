import sys
sys.path.append("..")

import booleantools as bt
import json
import ast
import multiprocessing
from itertools import permutations
from booleantools import BooleanFunction

def main():
    permset = list(permutations([0,1,2,3,4]))
    f_set = []
    with open("in/all", "r") as file:
        for line in file:
            f = BooleanFunction(ast.literal_eval(line), 5)
            f_set.extend(bt.perms_orbit_polynomial(permset, f))
    
    #Split f_set into 20 parts
    n = len(f_set) // 20
    packed = [f_set[start:start+n] for start in range(0, len(f_set), n)]

    queue = multiprocessing.Queue()
    
    #Start threads
    threads = [multiprocessing.Process(target=process_dataset, args=(dataset,queue)) for dataset in packed]
    
    [proc.start() for proc in threads]

    file_output(queue, len(threads))

def process_dataset(dataset, queue):
    for f in dataset:
        queue.put(bt._bin_to_dec(f.tableform))

    queue.put("DONE")

def file_output(queue, numThreads):
    with open("ruleNums.csv", "w") as file:
        should_stop = False
        done_count = 0
        while not should_stop:
            if not queue.empty():
                ruleNo = queue.get_nowait()
                if ruleNo == "DONE":
                    done_count += 1
                    if done_count == numThreads:
                        should_stop = True

                elif ruleNo is not None:
                    file.write(str(ruleNo) + ",")



main()
