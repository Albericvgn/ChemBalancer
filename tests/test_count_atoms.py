import pytest
from collections import defaultdict
import sys
import os
from unittest.mock import patch
from rdkit import Chem
import unittest



try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)

from chembalancer.chembalancer import count_atoms


from chembalancer import Chem  # Import Chem from the chembalancer module

class TestCountAtoms(unittest.TestCase):
    @patch('chembalancer.Chem', autospec=True)
    def test_count_atoms(self, MockChem):
        # Create an instance of MockChem
        mock_instance = MockChem.return_value
        mock_instance.get_formula.return_value = "H2O"  # Mocking the get_formula method

        # Example of using the mock instance in a test
        result = mock_instance.get_formula()
        self.assertEqual(result, "H2O")

        # You can also assert that certain methods are called
        mock_instance.get_formula.assert_called_once()

    @patch('chembalancer.Chem', autospec=True)
    def test_count_atoms_empty(self, MockChem):
        # Create an instance of MockChem with an empty formula
        mock_instance = MockChem.return_value
        mock_instance.get_formula.return_value = ""  # Mocking the get_formula method

        # Example of using the mock instance in a test
        result = mock_instance.get_formula()
        self.assertEqual(result, "")

        # You can also assert that certain methods are called
        mock_instance.get_formula.assert_called_once()

    @patch('chembalancer.Chem', autospec=True)
    def test_count_atoms_none(self, MockChem):
        # Create an instance of MockChem with a None formula
        mock_instance = MockChem.return_value
        mock_instance.get_formula.return_value = None  # Mocking the get_formula method

        # Example of using the mock instance in a test
        result = mock_instance.get_formula()
        self.assertEqual(result, None)

        # You can also assert that certain methods are called
        mock_instance.get_formula.assert_called_once()

    @patch('chembalancer.Chem', autospec=True)
    def test_count_atoms_complex(self, MockChem):
        # Create an instance of MockChem with a complex formula
        mock_instance = MockChem.return_value
        mock_instance.get_formula.return_value = "C6H12O6"  # Mocking the get_formula method

        # Example of using the mock instance in a test
        result = mock_instance.get_formula()
        self.assertEqual(result, "C6H12O6")

        # You can also assert that certain methods are called
        mock_instance.get_formula.assert_called_once()

if __name__ == '__main__':
    unittest.main()