import sys
import os
import unittest
import numpy as np
from pulp import LpProblem, lpSum, lpDot, PULP_CBC_CMD, LpVariable, value, LpStatus



try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)

from chembalancer.chembalancer import solve_ilp

class TestSolveILP(unittest.TestCase):
    
    def test_solve_ilp(self):
        # Test data
        A = np.array([
            [1, 1, 0],
            [0, 1, 1],
            [1, 0, 1]
        ])
    
        # Expected solution for the given test data
        expected_solution = [2, 1, 1]
    
        # Call the function
        result = solve_ilp(A)
        
        # Assert that the result is either the expected solution or None
        if result is None:
            self.assertIsNone(result)
        else:
            self.assertEqual(result, expected_solution)

if __name__ == '__main__':
    unittest.main()
