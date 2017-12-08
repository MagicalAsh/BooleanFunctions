import sys
sys.path.append("..")

import random
import booleantools as bt

def main():
    lst = [[0], [1], [2], [3], [4], [0,1]]
    rule = bt.BooleanFunction(lst, 5)
    back = [random.randint(0,1) for i in range(0,190)]
    back[80] = 1
    
    time = []
    ca = bt.CellularAutomata(rule, back)

    for i in range(0,300):
        time.append(ca.state)
        pretty_print(ca.state)
        ca.update()

    print("CA Center Col %1's: ", sum(get_col(time, 80))/len(back))

def pretty_print(state):
    for cell in state:
        char = "\033[7m \033[0m" if cell == 1 else " "
        print(char, sep="", end="")

    print()

def get_col(state, colNo):
    return [row[colNo] for row in state]

if __name__ == "__main__":
    main()
