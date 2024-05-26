import unittest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sys
import os

try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)

from chembalancer.chembalancer import get_molecular_formula



class TestChemFunctions(unittest.TestCase):
    def test_valid_smiles(self):
        # Test with a valid SMILES string for water
        self.assertEqual(get_molecular_formula('O'), 'H2O')

    def test_invalid_smiles(self):
        # Test with an invalid SMILES string
        self.assertEqual(get_molecular_formula('XYZ'), 'Invalid SMILES string')

    def test_another_valid_smiles(self):
        # Test with a valid SMILES string for benzene
        self.assertEqual(get_molecular_formula('C1=CC=CC=C1'), 'C6H6')

if __name__ == '__main__':
    unittest.main()
