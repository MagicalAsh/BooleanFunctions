import booleantools as bt

def main():
    r30 = bt.BooleanFunction([[0],[1],[2],[1,2]], 3)
    ca = bt.CellularAutomata(r30, [1,0,0,0,0])

    for i in range(0,3):
        print(ca.state)
        ca.update()

if __name__ == "__main__":
    main()
