
from src/chembalancer.chembalancer import hello_smiles


# Test the function
def test_hello_smiles():
    assert hello_smiles("C(=O)O") == "Hello, C(=O)O!", "Test failed: SMILES input"
