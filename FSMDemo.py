import booleantools
from booleantools import BooleanFunction
from booleantools import FiniteStateMachine

def main():
    r34 = BooleanFunction([[0,2],[1,2],[0,1],[1]], 3)
    fsm = FiniteStateMachine(r34)
    for row in fsm.adjMatr:
        print(row)

main()
