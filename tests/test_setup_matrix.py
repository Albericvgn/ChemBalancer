import numpy as np
import sys
import os
import pytest

try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)


from chembalancer.chembalancer import setup_matrix

def test_setup_matrix():
    elements = ['H', 'C', 'O']
    counts = [
        {'H': 2, 'C': 1},
        {'H': 2, 'O': 1},
        {'C': 1, 'O': 1},
    ]
    
    expected_matrix = np.array([
        [2, 1, 0],
        [2, 0, 1],
        [0, 1, 1]
    ])
    
    matrix = setup_matrix(elements, counts)
    
    assert isinstance(matrix, np.ndarray), "Output should be a numpy array"
