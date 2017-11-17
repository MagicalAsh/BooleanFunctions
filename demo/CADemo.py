import sys
sys.path.append("..")

import booleantools as bt

def main():
    lst = [[0],[1],[2],[3], [2,3]]
    rule = bt.BooleanFunction(lst, 5)
    back = [0 for i in range(0,150)]
    back[75] = 1
    
    ca = bt.CellularAutomata(rule, back)

    for i in range(0,300):
        pretty_print(ca.state)
        ca.update()

def pretty_print(state):
    for cell in state:
        char = "\033[7m \033[0m" if cell == 1 else " "
        print(char, sep="", end="")

    print()

if __name__ == "__main__":
    main()
