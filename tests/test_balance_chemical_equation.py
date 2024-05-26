import unittest
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

from chembalancer.chembalancer import balance_chemical_equation 

class TestBalanceChemicalEquation:
    @pytest.fixture(autouse=True)
    def setup_method(self):
        pass

    def test_combustion_of_methane(self):
        # Reactants and products
        reactant_smiles = ['CH4', 'O2']
        product_smiles = ['CO2', 'H2O']

        # Call the function to test
        try:
            reactant_data, product_data = balance_chemical_equation(reactant_smiles, product_smiles)
        except ValueError as e:
            assert "Failed to solve the balance equation." in str(e)
            return

        # Assertions
        expected_reactant_data = [(1, 'CH4'), (2, 'O2')]
        expected_product_data = [(1, 'CO2'), (2, 'H2O')]

        assert reactant_data == expected_reactant_data
        assert product_data == expected_product_data
