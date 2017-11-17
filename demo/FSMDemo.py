import sys
sys.path.append("..")

import json
import booleantools
from booleantools import BooleanFunction
from booleantools import FiniteStateMachine
import copy

def main():
    with open("in/2res_classes.json") as fileIn:
        f_list = [BooleanFunction(f, 5) for f in json.load(fileIn)["classes"]]
    
    for f in f_list:
        fsm = FiniteStateMachine(f)
        for line in fsm.adjMatr:
            print(*line, sep="\t")
        words = search(fsm, 0, 5)
        legitWords = set([])
        for word in getSubVectors(5, words):
            tupledWord = tuple(word)
            legitWords.add(tupledWord)

        print("Words = ", legitWords)
        print("\n\n")

def getSubVectors(length, listStuff):
    for i in listStuff:
        if len(i) != length:
            yield from getSubVectors(length, i)
        elif isinstance(i[0], int):
            yield i

def search(fsm, start, numMore, pathTaken=[]):
    if numMore == 0: 
        return pathTaken
    
    out = []
    for i in range(0, len(fsm.adjMatr)):
        adj = [(x,weight) for x,weight in enumerate(fsm.adjMatr[i]) if weight is not None]
        for i in adj:
            newRow = fsm.adjMatr[i[0]]
            newList = copy.deepcopy(pathTaken) + [i[1]]
            out.extend([search(fsm, newRow, numMore - 1, newList)])
        
    return out 


main()
