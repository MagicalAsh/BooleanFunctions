import booleantools as bt


r30 = bt.BooleanFunction([[0],[1],[2],[1,2]], 3)

tca = bt.TargetedCellularAutomata(r30, [1,1,0,1,0], [0,0,1,1,0], *bt.TargetTransitions.standard_revs()) 

tca.pretty_print()
print("\n\n")

[tca.update() for i in range(1,10)]
tca.pretty_print()

for row in tca.time:
    print(row)
