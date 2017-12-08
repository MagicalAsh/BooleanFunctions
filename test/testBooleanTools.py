import sys
sys.path.append("..")

import booleantools as bt
import unittest

r30 = bt.BooleanFunction([[0],[1],[2],[1,2]], 3)
r90 = bt.BooleanFunction([[0],[2]], 3)

class BTUnitTest(unittest.TestCase):
    
    def test_hamming_wt(self):
        self.assertEqual(r30.hamming_weight(), 4)

    def test_ham_dist(self):
        self.assertEqual(r30.hamming_distance(r90), 2)
        self.assertEqual(r30.hamming_distance([r90, r30]), [2, 0])
        
    def test_nonlinearity(self):
        self.assertEqual(r30.nonlinearity(), 2)
        self.assertEqual(r90.nonlinearity(), 0)

    def test_eval(self):
        r30p1 = bt.BooleanFunction([[0],[1],[2],[1,2],[]], 3)
        self.assertEqual(r30(1,0,0), 1)
        self.assertEqual(r30(0,0,0), 0)
        self.assertEqual(r30p1(0,0,0), 1)

        self.assertEqual(r90(1,0,1), 0)
        self.assertEqual(r90(1,1,0), 1)

    def test_perm(self):
        self.assertEqual(r30.apply_permutation([0,1,2]), r30)
        self.assertNotEqual(r30.apply_permutation([1,0,2]), r30)
        
    def test_balanced(self):
        self.assertTrue(r30.is_balanced())

    def test_correlation_immune(self):
        self.assertFalse(r30.is_correlation_immune())
        self.assertTrue(r90.is_correlation_immune())
    
    def test_linear_structs(self):
        self.assertEqual(r30.linear_structures(), set([0]))
        self.assertEqual(r90.linear_structures(), set([0,2]))
    
    def test_str(self):
        self.assertEqual(repr(r30), "BooleanFunction([[0], [1], [2], [1, 2]], 3)")
        self.assertEqual(str(r30), "x_{0} \\oplus x_{1} \\oplus x_{2} \\oplus x_{1}x_{2}")
        self.assertEqual(r30.tex_str(True), "$x_{0} \\oplus x_{1} \\oplus x_{2} \\oplus x_{1}x_{2}$")
        
    def test_add(self):
        rule = bt.BooleanFunction([[1],[1,2]], 3)
        self.assertEqual(r30 + r90, rule)

    def test_hash(self):
        self.assertEqual(hash(r30), 30)
        self.assertEqual(hash(r90), 90)

    def test_k_resilient(self):
        self.assertFalse(r30.is_k_resilient())
        self.assertTrue(r90.is_k_resilient())
        self.assertFalse(r90.is_k_resilient(k=2))

    def test_FSM(self):
        fsm = bt.FiniteStateMachine(r30)
        self.assertEqual(fsm.adjMatr[0][0], 0)
        self.assertEqual(fsm.adjMatr[3][0], None)

    def test_CA(self):
        ca = bt.CellularAutomata(r30, [0,1,1,0,0])
        ca.update()
        self.assertEqual(ca.state[1], 1)
        self.assertEqual(ca.get_column(1), [1,1])

if __name__ == "__main__":
    unittest.main()
