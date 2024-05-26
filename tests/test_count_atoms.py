import pytest
from collections import defaultdict

try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)

from chembalancer.chembalancer import count_atoms  

# Mocking RDKit's Chem module for testing purposes
class MockAtom:
    def __init__(self, symbol):
        self.symbol = symbol
    
    def GetSymbol(self):
        return self.symbol

class MockMol:
    def __init__(self, atoms):
        self.atoms = atoms
    
    def GetAtoms(self):
        return self.atoms

def mock_mol_from_smiles(smiles):
    atom_counts = defaultdict(int)
    atoms = [MockAtom(symbol) for symbol in smiles]
    return MockMol(atoms)

# Test case for counting atoms
def test_count_atoms():
    smiles = "CC(=O)O"
    expected_result = {'C': 2, 'O': 2}
    
    # Mocking Chem.MolFromSmiles with a custom function
    count_atoms._mock_mol_from_smiles = mock_mol_from_smiles
    
    result = count_atoms(smiles)
    assert result == expected_result

# Test case for an empty SMILES string
def test_count_atoms_empty():
    smiles = ""
    expected_result = {}
    
    # Mocking Chem.MolFromSmiles with a custom function
    count_atoms._mock_mol_from_smiles = mock_mol_from_smiles
    
    result = count_atoms(smiles)
    assert result == expected_result

# Test case for None input
def test_count_atoms_none():
    smiles = None
    expected_result = {}
    
    # Mocking Chem.MolFromSmiles with a custom function
    count_atoms._mock_mol_from_smiles = mock_mol_from_smiles
    
    result = count_atoms(smiles)
    assert result == expected_result

# Test case for a complex SMILES string
def test_count_atoms_complex():
    smiles = "C1=CC=CC=C1C(=O)O"
    expected_result = {'C': 7, 'H': 6, 'O': 2}
    
    # Mocking Chem.MolFromSmiles with a custom function
    count_atoms._mock_mol_from_smiles = mock_mol_from_smiles
    
    result = count_atoms(smiles)
    assert result == expected_result

if __name__ == "__main__":
    pytest.main()
