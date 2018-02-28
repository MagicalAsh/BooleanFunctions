import sys
sys.path.append("..")
import ast
import json
import random
import booleantools as bt
from booleantools import BooleanFunction, CellularAutomata

def main():
    
    with open("in/grouped/card120", "r") as filein:
        f_list = [BooleanFunction(f, 5) for f in [ast.literal_eval(line) for line in filein]]
    
    for ca in f_list:
        if ca.linear_structures() != set():
            time = evalCA(ca)

    
def evalCA(rule):
    ca = getCA(rule)
    print("Testing CA ", rule)
    print("Column balanced: ", colBalanced(ca))
    print("Order Balanced: ", normalize(orderBalanced(ca)))
    print("\n")


def getCA(rule, iterations=300):
    time = []
    state = [random.randint(0,1) for i in range(150)]
    ca = CellularAutomata(rule, state)

    for i in range(0, iterations):
        time.append(ca.state)
        ca.update()

    return time
    
def colBalanced(ca):
    col = get_col(ca, 75)
    return sum(col)/len(col)

def orderBalanced(ca):
    dic = {0: {0:0, 1:0}, 1:{0:0, 1:0}}
    col = get_col(ca, 75)
    for i in range(len(col) - 1):
        here = col[i]
        nex = col[i+1]
        dic[here][nex] += 1

    return dic

def normalize(dic):
    newDic = {}
    for key in dic:
        newDic[key] = {key2: dic[key][key2]/len(dic[key].values()) for key2 in dic[key]}

    return newDic

def pretty_print(state):
    for cell in state:
        char = "\033[7m \033[0m" if cell == 1 else " "
        print(char, sep="", end="")

    print()

def get_col(state, colNo):
    return [row[colNo] for row in state]

if __name__ == "__main__":
    main()
